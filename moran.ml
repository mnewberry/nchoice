
(* A Pop (population) is the state datastructure, which counts types. *)

module BS = BatSet

module type P = sig
  type 'a pop
  val empty : 'a pop
  val size : 'a pop -> int
  val max_pop : 'a pop -> int
  val max_type : 'a pop -> 'a
  val add: int -> 'a -> 'a pop -> 'a pop
  val del: int -> 'a -> 'a pop -> 'a pop
  val delf : (unit -> unit) -> int -> 'a -> 'a pop -> 'a pop
  val to_list : 'a pop -> (int * int * 'a) list
  val of_list : ('a * int) list -> 'a pop
  val to_set_f : ((int * int * 'a) -> 'b) -> 'a pop -> 'b BS.t
  val to_assoc : 'a pop -> ('a * int) list
  val lookup : 'a -> 'a pop -> int
  val to_str : ('a -> string) -> 'a pop -> string
  val nth_indiv : int -> 'a pop -> 'a
  val choose : 'a pop -> 'a
  val fold : (int -> 'a -> 'b -> 'b) -> 'b -> 'a pop -> 'b
  val fold_indiv : ('a -> 'b -> 'b) -> 'b -> 'a pop -> 'b
end

module Pop : P = struct
  type 'a pop = Empty | Node of (int * int * 'a * 'a pop)

  let empty = Empty

  let tot_pop = function
      Empty -> 0
    | Node (tpop, _, _, _) -> tpop

  let size = tot_pop

  let max_pop = function
      Empty -> 0
    | Node (_, npop, _, _) -> npop

  let max_type = function
      Empty -> raise Not_found
    | Node (_, _, ty, _) -> ty

  let expect_node = function
      Empty -> raise Not_found
    | Node data -> data

  let rec add n ty pop = if n = 0 then pop else match pop with 
      Empty -> Node (n, n, ty, Empty)
    | Node (ap, bp, bty, rest) ->
       if ty = bty then Node (ap + n, bp + n, bty, rest) else
       let (nrap, nrbp, nrbty, nrrest) as nr = expect_node (add n ty rest) in
       if nrbp > bp
       then Node (nrap + bp, nrbp, nrbty, 
              Node (ap - (nrbp - n), bp, bty, nrrest))
       else Node (ap + n, bp, bty, Node nr)

  let rec insert (p, ty) = function
      Empty -> Node (p, p, ty, Empty)
    | Node (ap, bp, bty, rest) ->
        if p >= bp then Node (p + ap, p, ty, Node (ap, bp, bty, rest))
        else Node (ap + p, bp, bty, insert (p, ty) rest)

  let rec delf f n ty = function
      Empty -> raise Not_found
    | Node (ap, bp, bty, rest) ->
       if ty <> bty then Node (ap - n, bp, bty, delf f n ty rest) else 
       if bp - n < 0 then raise Not_found else 
       if bp - n = 0 then (f ; rest) else insert (bp - n, bty) rest 

  let del n ty = delf (fun () -> ()) n ty

  let of_list l = Mu.fold (fun (ty, n) pop -> add n ty pop) empty l

  let rec to_list = function
      Empty -> []
    | Node (ap, bp, bty, rest) -> (ap, bp, bty) :: to_list rest

  let rec to_set_f f = function
      Empty -> BS.empty
    | Node (ap, bp, bty, rest) -> BS.add (f (ap, bp, bty)) (to_set_f f rest)

  let rec to_assoc = function
      Empty -> []
    | Node (ap, bp, bty, rest) -> (bty, bp) :: to_assoc rest

  let rec lookup ty = function
      Empty -> raise Not_found
    | Node (_, bp, bty, r) -> if bty = ty then bp else lookup ty r

  let to_str f state =
    Printf.sprintf "[%s]" (String.concat "; "
      (Mu.map (fun (an, bn, ty) -> Printf.sprintf "%dx%s" bn (f ty))
              (to_list state)))

  (* one-based, the ordering is obscure *)
  let rec nth_indiv n = function
      Empty -> raise Not_found
    | Node (ap, bp, bty, rest) ->
       if n <= 0 then raise Not_found else 
       if n <= bp then bty else nth_indiv (n - bp) rest

  let choose p = nth_indiv (Rand.int (size p) + 1) p

  let rec fold kons knil = function
      Empty -> knil
    | Node (ap, bp, bty, rest) -> fold kons (kons bp bty knil) rest

  (* call kons once for each individual, rather than each genotype *)
  let fold_indiv kons knil lyst =
    let nk n genotype kn = Mu.rec_n n (kons genotype) kn in
    fold nk knil lyst
end

(* fitness is a Meas (measure) on genotypes *)
module Meas = struct
  type 'a meas = Empty | Node of (float * float * 'a * 'a meas)

  let empty = Empty

  let total = function Empty -> 0. | Node (a, _, _, _) -> a
  let expect_node = function Empty -> raise Not_found | Node data -> data

  let rec fadd x ty meas = if x = 0. then meas else match meas with
      Empty -> Node (x, x, ty, Empty)
    | Node (ap, bp, bty, rest) ->
       if ty = bty then Node (ap +. x, bp +. x, bty, rest) else
       let (nrap, nrbp, nrbty, nrrest) as nr = expect_node (fadd x ty rest) in
       if nrbp > bp
       then Node (nrap +. bp, nrbp, nrbty, 
              Node (ap -. (nrbp -. x), bp, bty, nrrest))
       else Node (ap +. x, bp, bty, Node nr)

  let rec insert (p, ty) = function
      Empty -> Node (p, p, ty, Empty)
    | Node (ap, bp, bty, rest) ->
        if p >= bp then Node (p +. ap, p, ty, Node (ap, bp, bty, rest))
        else Node (ap +. p, bp, bty, insert (p, ty) rest)

  let rec elt x = function
      Empty -> raise Not_found
    | Node (ap, bp, bty, rest) ->
       if x < 0. then raise Not_found else
       if x <= bp then bty else elt (x -. bp) rest

  let choose m = elt (Rand.float (total m)) m

  let rec fold kons knil = function
      Empty -> knil
    | Node (ap, bp, bty, rest) -> fold kons (kons bp bty knil) rest

  let merge inp meas = fold fadd meas inp

  (* multinomially sample the population *)
  let sample_pop pop_size meas =
    let count bp bty (t_k, p_k, n) = (bty :: t_k, bp :: p_k, n + 1) in
    let (t_k, p_k, len) = fold count ([],[],0) meas in
    let n_k = Rand.multinom pop_size p_k in
    Mu.fold2 Pop.add Pop.empty n_k t_k

  let rec to_list = function
      Empty -> []
    | Node (ap, bp, bty, rest) -> (ap, bp, bty) :: to_list rest

  let of_list l = Mu.fold (fun (ty, n) pop -> fadd n ty pop) empty l

  let rec to_alist meas = fold (fun bp bty knil -> (bty, bp) :: knil) [] meas

  let rec measure ty = function
      Empty -> 0.
    | Node (ap, bp, bty, rest) -> if bty = ty then bp else measure ty rest

  let prob ty = function
      Empty -> 0.
    | (Node (ap, bp, bty, rest) as nd) -> measure ty nd /. ap

  let to_str f meas =
    Printf.sprintf "[%s]" (String.concat "; "
      (Mu.map (fun (an, bn, ty) -> Printf.sprintf "%fx%s" bn (f ty))
              (to_list meas)))

  let rec rescale factor = function
      Empty -> Empty
    | Node (ap, bp, bty, rest) -> 
        Node (ap *. factor, bp *. factor, bty, rescale factor rest)

  let normalize = function Empty -> Empty 
    | Node (ap,_,_,_) as node -> rescale (1. /. ap) node
 
end
