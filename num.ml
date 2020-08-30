

(**** VECTORs *********)
type 'a vector = Vec of 'a list

(**** constructors ****)
let nullvec n =
  let rec inner (n, acc) =
    if n = 0 then acc else inner (n-1, 0. :: acc)
in
  Vec (inner (n, []))

let basevec n k =
  let rec inner (i, acc) =
    if i = n+1 then Vec acc
  else inner (i+1, (if i = n-k+1 then 1. else 0.)::acc)
in
  inner (1, [])

let random_vec n =
  Random.self_init ();
  let rec inner n acc =
    if n = 0 then acc
  else inner (n-1) ((Random.float 1.)::acc)
in
  Vec (inner n [])


(**** functions *******)
let tl_vec (Vec v) = Vec (List.tl v)
let hd_vec (Vec v) = Vec (List.hd v)

let dot_vec (Vec a) (Vec b) =
  let rec inner acc = function
    | [], [] -> acc
    | x::xs, y::ys -> inner (acc +. x *. y) (xs, ys)
    | _  -> failwith "Length mismatch"
in
  inner 0. (a, b)

let plus_vec (Vec a) (Vec b) =
  Vec (List.map2 (+.) a b)

let scalar_mult_vec lambda (Vec v) = Vec (List.map (( *.) lambda) v)

(****** MATRICes *******
 as list of row vecs, list of column vecs *)
type 'a matrix = Mat of ('a vector vector * 'a vector vector)

(**** constructors ****)
let matrix_from_rows rows =
  let rec inner (acc, rem) = try ( inner (
      Vec (List.fold_right List.cons (List.map List.hd rem) []) :: acc,
      List.map List.tl rem
    ) ) with e -> Vec (List.rev acc)
in
  Mat (Vec (List.map (fun r -> Vec r) rows), inner ([], rows))

let identity n =
  let rec inner acc k =
    if k = n+1 then Vec acc
  else inner (basevec n (n-k+1) :: acc) (k+1)
in
  let bases = inner [] 1
in
  Mat (bases, bases)

let nullmatrix n =
  let rec inner acc k =
    if k = 0 then Vec acc
  else inner (nullvec n :: acc) (k-1)
in
  let nullvecs = inner [] n
in
  Mat (nullvecs, nullvecs)

let random_mat n =
  let rec inner k acc =
    if k = 0 then acc
  else
    let (Vec rv) = random_vec n
  in
    inner (k-1) (rv::acc)
in
  inner n [] |> matrix_from_rows


(**** functions *******)
let matrix_to_lists (Mat (Vec rowvecs, _)) = List.map (fun (Vec r) -> r) rowvecs

let transpose (Mat (rows, cols)) = Mat (cols, rows)

let dot_mat_vec (Mat (Vec rows, _)) x = Vec (List.map (dot_vec x) rows)

let plus_mat (Mat (Vec ar, Vec ac)) (Mat (Vec br, Vec bc)) =
  Mat (
    Vec (List.map2 plus_vec ar br),
    Vec (List.map2 plus_vec ac bc)
  )

let dot_mat (Mat (ar, ac)) (Mat (br, bc)) =
  let Vec arl, Vec bcl = ar, bc
in
  Mat (
    Vec (List.map (dot_mat_vec (Mat (bc, Vec []))) arl),
    Vec (List.map (dot_mat_vec (Mat (ar, Vec []))) bcl)
  )

let outer_prod (Vec u) (Vec v) =
  Mat (
    Vec (List.map (fun ui -> scalar_mult_vec ui (Vec v)) u),
    Vec (List.map (fun vi -> scalar_mult_vec vi (Vec u)) v)
  )

let rec list_drop n ls = if n = 0 then ls else match ls with
  | _::xs -> list_drop (n-1) xs
  | [] -> failwith "Drop more than length"

let list_take n ls =
  let rec inner (n, acc) = if n = 0 then fun _ -> acc else function
    | [] -> failwith "Take more than length"
    | x::xs -> inner (n-1, x::acc) xs
in
  inner (n, []) ls


let prepend_zeros n (Vec v) =
  let rec inner (n, v) =
    if n = 0 then v else inner (n-1, 0.::v)
in
  Vec (inner (n, v))






(* CHOLESKY DECOMPOSITION *)

(* O(1) *)
let cholesky_get_a_b_B (Mat (Vec rows, Vec cols)) =
  let (Vec first) = List.hd cols
in
  let (a, b) = match first with
    | [] -> failwith "No more row to cholesky_get"
    | a::b -> (a, Vec b)
in
  let bmat =
    let brows = Vec (List.map tl_vec (List.tl rows))
  in
    let bcols = Vec (List.map tl_vec (List.tl cols))
  in
    Mat (brows, bcols)
in
  ((a, b), bmat)


(* O(n^2) *)
let cholesky_get_L n i (a, Vec b) =
  let root_a = assert (a > 0.); sqrt a
in
  let Mat (Vec idrows, Vec idcols) = identity (n-i)
in
  let idcols = List.map (prepend_zeros i) idcols
in
  let (Vec scaled_b) = scalar_mult_vec (1. /. root_a) (Vec b)
in
  let rows = List.map2 (fun bi -> fun (Vec idrow) -> Vec (bi::idrow)) scaled_b idrows
in
  let rows = List.map (prepend_zeros (i-1)) rows
in
  let abcol = prepend_zeros (i-1) (Vec (root_a :: scaled_b))
in
  let cols =
    let rec loop (rem, acc) =
      if rem = 0 then Vec acc
    else loop (rem-1, (basevec n rem)::acc)
  in
    loop (i-1, abcol::idcols)
in
  let rows = (scalar_mult_vec root_a (basevec n i))::rows
in
  let rows =
    let rec loop (rem, acc) =
      if rem = 0 then Vec acc
    else loop (rem-1, (basevec n rem)::acc)
  in
    loop (i-1, rows)
in
  (* let matlist = (Mat (rows, cols)) |> transpose |> matrix_to_lists in ignore (List.map (fun lis -> print_string (String.concat " " (List.map string_of_float lis)); print_newline ()) matlist); print_newline (); *)
  Mat (rows, cols)

(* wrapper *)
let cholesky_get_next_a_L n i amat =
  let ((a, b), bmat) = cholesky_get_a_b_B amat
in
  let scaled_bb = outer_prod (scalar_mult_vec (-1. /. a) b) b
in
  let next_amat = plus_mat bmat scaled_bb
in
  let lmat = cholesky_get_L n i (a, b)
in
  (next_amat, lmat)

(* O(n^3) *)
let cholesky_decomp mat =
  let (Mat (Vec rows, Vec cols)) = mat
in
  let n = let l = List.length rows in assert (l = List.length cols); l
in
  let rec loop (i, l_acc, amat) =
    if i = n+1 then l_acc
  else
    let (next_amat, lmat) = cholesky_get_next_a_L n i amat
  in
    loop (
      i+1,
      dot_mat l_acc lmat,
      next_amat
    )
in
  loop (1, identity n, mat)


(* Forward substitution 
 i.e. solution to Lx = b *)

let forward_substitute (Mat (Vec lrows, _)) (Vec b) =
  let rec partial_dot acc = function
    | l::ls, x::xs -> partial_dot (acc +. l *. x) (ls, xs)
    | l::ls, [] -> (l, acc)
    | _, _ -> failwith "Too long partial dot"
in
  let rec loop xacc = function
    | [], [] -> Vec xacc
    | b::bs, (Vec lrow)::lrows ->
        let (l, prod) = partial_dot 0. (lrow, xacc)
      in
        loop (xacc @ [(b -. prod) /. l]) (bs, lrows)
    | _, _ -> failwith "Length mismatch"
in
  loop [] (b, lrows)


(* Backward substitution
 i.e. solution to Ux = b *)
let backward_substitute (Mat (Vec urows, _)) (Vec b) =
  let n = let n = List.length b in assert ((List.length urows) = n); n
in
  let rec loop i xacc = function
  | [], [] -> Vec xacc
  | b::bs, (Vec urow)::urows ->
    let (u, utail) = match (list_drop (n-i) urow) with
      | u::utail -> (u, utail)
      | [] -> failwith "u row too short"
  in
    let prod = dot_vec (Vec utail) (Vec xacc)
  in
    loop (i+1) (((b -. prod) /. u)::xacc) (bs, urows)
  | _, _ -> failwith "Length mismatch"
in
  loop 1 [] (List.rev b, List.rev urows)


(* Solution to linear system Mx = b
   with precomputing the cholesky decomp so that it's faster for multiple b's *)
let solve_linear_system mat =
  let l = cholesky_decomp mat
in
  let u = transpose l
in
  fun b -> forward_substitute l b |> backward_substitute u


(* Same but without precomputing the decomposition to allow O(1) currying *)
let solve_linear_system2 mat b =
  let l = cholesky_decomp mat
in
  let u = transpose l
in
  forward_substitute l b |> backward_substitute u




(***** STREAMs ********)
type 'a stream = Cons of 'a * (unit -> 'a stream)

(**** constructors ****)
let rec const k = Cons (k, fun () -> const k)

(**** functions *******)
let hd (Cons (h, _)) = h
let tl (Cons (_, t)) = t

let rec take n (Cons (v, f)) =
  if n = 0 then [] else take (n-1) (f ()) |> List.cons v


let limit seq =
  let rec inner prev (Cons (v, f)) =
    if v -. prev |> Float.abs < 1e-10 then v
  else inner v (f ())
in
  let (Cons (v, f)) = seq
in
  inner v (f ())


let limit_vec seq =
  let rec inner prev (Cons (Vec v, f)) =
    if (List.for_all2 (fun vi -> fun previ -> (vi -. previ) < 1e-10) v prev) then Vec v
  else inner v (f ())
in
  let (Cons (Vec v, f)) = seq
in
  inner v (f ())




(* FUNCTIONS R^n -> R *)

(* derivatives *)
(* R -> R *)
let d_1dim f x =
  let rec seq delta = Cons (
    (f (x +. delta) -. f (x -. delta)) /. (2. *. delta),
    fun () -> seq (delta /. 2.)
  )
in
  seq 0.1 |> limit

(* partial derivative with respect to ith component, 1 indexed *)
(* each evaluation of substituting is O(i)! *)
let fix_all_but i (Vec x0) =
  let rec inw k acc = function
    | x0::x0s ->
      if k = 0 then ((fun x -> outw (x::x0s) acc), x0) else inw (k-1) (x0::acc) x0s
    | _ -> failwith "Length mismatch"
  and outw x0s = function
    | ac::acs -> outw (ac::x0s) acs
    | _ -> Vec x0s
in
  inw (i-1) [] x0

let partial i f (Vec x0) =
  let (f_1dim, x0i) = fix_all_but i (Vec x0)
in
  d_1dim (fun x -> f_1dim x |> f) x0i


(* total derivative *)
let d n f x =
  let rec loop i acc =
    if i = 0 then Vec acc
  else loop (i-1) ((partial i f x)::acc)
in
  loop n []

(* Hessian *)
let d_1_to_n n f x =
  let rec seq delta = Cons (
    (plus_vec (f (x +. delta)) (scalar_mult_vec (-1.) (f (x -. delta)))) |> scalar_mult_vec (1. /. (2. *. delta)),
    fun () -> seq (delta /. 2.)
  )
in
  seq 0.1 |> limit_vec


let partial_vec n i f (Vec x0) =
  let (f_1_to_n, x0i) = fix_all_but i (Vec x0)
in
  d_1_to_n n (fun x -> f_1_to_n x |> f) x0i

let d2 n f x =
  let rec loop i acc =
    if i = 0 then Vec acc
  else loop (i-1) ((partial_vec n i (d n f) x)::acc)
in
  let dij = loop n []
in
  Mat (dij, dij)



(* GRADIENT DESCENT *)
let graddesc_minimise n f x0 delta =
  let rec seq x = Cons (
    x,
    fun () -> d n f x |> scalar_mult_vec (-1. *. delta) |> plus_vec x |> seq
  )
in
  seq x0 |> limit_vec

(* NEWTON'S METHOD *)
let newton_minimise n f x0 =
  let rec seq x = Cons (
    x,
    fun () ->
      let h = solve_linear_system (d2 n f x) (scalar_mult_vec (-1.) (d n f x))
    in
      seq (plus_vec x h)
  )
in
  seq x0 |> limit_vec








let testfun x = match x with
  | Vec (x::y::z::[]) -> x**3.*.z +. y*.z**2.
  | _ ->  failwith "wrong test"

let test = d 3 testfun (Vec [1.; 2.; 3.])
let test2 = d2 3 testfun (Vec [1.; 2.; 3.])



let testfun2 x = match x with
  | Vec (x1::x2::[]) -> 4.*.x1*.x1 +. x2*.x2 -. 2.*.x1*.x2
  | _ -> failwith "wrong test"



let lsol = solve_linear_system (matrix_from_rows [
  [4.; 12.; -16.];
  [12.; 37.; -43.];
  [-16.; -43.; 98.]
]) (
  Vec [2.; 5.; -8.]
)

let lsol2 = solve_linear_system2 (matrix_from_rows [
  [4.; 12.; -16.];
  [12.; 37.; -43.];
  [-16.; -43.; 98.]
]) (
  Vec [2.; 5.; -8.]
)


let usol = backward_substitute (matrix_from_rows [
  [2.; -1.; 0.];
  [0.; 3.; -2.];
  [0.; 0.; 1.];
]) (Vec [-1.; 7.; 1.])



let a = matrix_from_rows [[3.; 2.]; [-1.; 4.]];;
let b = matrix_from_rows [[-2.; 1.]; [3.; 3.]];;

let temp = cholesky_decomp (matrix_from_rows [
  [4.; 12.; -16.];
  [12.; 37.; -43.];
  [-16.; -43.; 98.]
])







