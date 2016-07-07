let fl = float_of_int
let trunc = truncate
type genotype = G_aa | G_aA | G_AA

(* select the G_??th value from a (xaa, xaA, xAA) vector *)
let popct (paa, paA, pAA) genotype = match genotype with 
    G_aa -> paa
  | G_aA -> paA
  | G_AA -> pAA

let poppred (paa, paA, pAA) genotype = match genotype with 
    G_aa -> (max 0 (paa - 1), paA, pAA)
  | G_aA -> (paa, max 0 (paA - 1), pAA)
  | G_AA -> (paa, paA, max 0 (pAA - 1))

let genotypes = (G_aa,G_aA,G_AA)
let genotypes_list = [G_aa;G_aA;G_AA]

(* Vector shmarithmetic *)
let vzero = (0., 0., 0.)
let ( ++ ) (a, b, c) (aa, bb, cc) = (a + aa, b + bb, c + cc)
let ( ++. ) (a, b, c) (aa, bb, cc) = (a +. aa, b +. bb, c +. cc)
let vmap f (a, b, c) = (f a, f b, f c)
let ( **. ) aa v = vmap ( ( *. ) aa) v (* scalar float multiplication *)
let vinv = vmap ((-) 0)
let ( -- ) a b = a ++ vinv b
let vpr (a, b, c) = Printf.sprintf "(%d,%d,%d)" a b c
let vprf (a, b, c) = Printf.sprintf "(%f,%f,%f)" a b c
let vsum (a, b, c) = a + b + c
let vmin (a, b, c) = Mu.fold min a [b;c]
let vsumf (a, b, c) = a +. b +. c
let fdiv x n = fl x /. fl n
let vftimes aa v = vmap (fun a -> aa *. fl a) v

let countarray nn = let arr = Array.make nn 0 in
  Array.iteri (fun i _ -> arr.(i) <- i) arr ; arr

(* Multinomially sample a float vector and return an int vector *)
let multinom3 nn (paa, paA, pAA) = 
  let ss = Gsl.Randist.multinomial Rand.rng nn [|paa;paA;pAA|] in
  (ss.(0), ss.(1), ss.(2))

(* Mendelian offspring distribution vectors *)
let mendel a b = match (a,b) with
      G_aa, G_aa -> (1., 0., 0.)
    | G_aa, G_aA -> (0.5, 0.5, 0.)
    | G_aa, G_AA -> (0., 1., 0.)
    | G_aA, G_aa -> (0.5, 0.5, 0.)
    | G_aA, G_aA -> (0.25, 0.5, 0.25)
    | G_aA, G_AA -> (0., 0.5, 0.5)
    | G_AA, G_aa -> (0., 1., 0.)
    | G_AA, G_aA -> (0., 0.5, 0.5)
    | G_AA, G_AA -> (0., 0., 1.)

(* mch nn n s f pop: one stochastic update of the mate choice model *)
let mch_ multinom3 nn n (waa, waA, wAA) f pop =
  let (p_aa, p_aA, p_AA) = vmap (fun x -> fdiv x nn) pop in
  let ( - ) = ( -. ) and ( + ) = ( +. ) and ( * ) = ( *. )
    and i = 1. and h = 0.5 and q = 0.25 and n = fl n
    and nmi = fl n -. 1. in
  let pp_aa = p_aa * (i - (i - p_aa)**n)
    + h * p_aa * p_aA * (i - p_aA) ** nmi
    + h * p_aa * p_aA * (i - p_aa) ** nmi
    + q * (i - (i - p_aA) ** n) * p_aA
  and pp_aA =
      h * p_aA * (i - p_aA) ** nmi * (p_aa + p_AA)
    + h * p_aA * (p_aa * (i - p_aa) ** nmi + p_AA * (i - p_AA) ** nmi)
    + p_aa * p_AA * ((i - p_AA) ** nmi + (i - p_aa) ** nmi)
    + h * p_aA * (i - (i - p_aA) ** n)
  and pp_AA = p_AA * (i - (i - p_AA)**n)
    + h * p_AA * p_aA * (i - p_aA) ** nmi
    + h * p_AA * p_aA * (i - p_AA) ** nmi
    + q * (i - (i - p_aA) ** n) * p_aA in
  let interp f a b = f *. a +. (1. -. f) *. b in
  let ec_aa, ec_aA, ec_AA = (
    waa *. (interp f pp_aa p_aa),
    waA *. (interp f pp_aA p_aA),
    wAA *. (interp f pp_AA p_AA)) in
  let norm x = x /. (ec_aa +. ec_aA +. ec_AA) in
  multinom3 nn (norm ec_aa, norm ec_aA, norm ec_AA)

let mch = mch_ multinom3
let mch_det =
  let round (xaa, xaA, xAA) = 
    let rm f = mod_float f 1.0 and fl = trunc and cl f = trunc (ceil f) in
    if rm xaA==0. && rm xaA==0. && rm xAA==0. then (fl xaa, fl xaA, fl xAA)
    else if rm xaA >= rm xaa && rm xaA >= rm xAA then (fl xaa, cl xaA, fl xAA)
    else if rm xaa > rm xaA && rm xaa > rm xAA then (cl xaa, fl xaA, fl xAA)
    else (fl xaa, fl xaA, cl xAA) in
  let multinom3 nn freqs = round (fl nn **. freqs) in
  mch_ multinom3

let sterr n q = 1.96 *. sqrt (q *. (1. -. q) /. fl n)

(* Execute copies of f in nprocs parallel processes and pkons the results
   of each onto pknil.  Once the result of pkons is finishedp, continue
   pkonsing results, but don't spawn any new batches.  f must return an
   integer. *)
let unix_proc_fold nprocs f pkons pknil finishedp = 
  let procs = Array.make nprocs 
    (None : (int * in_channel * Unix.file_descr * Unix.file_descr) option) in

  let open_proc () = 
    let inde, outde = Unix.pipe () in match Unix.fork () with
        0 -> let result = f () in 
          (output_binary_int (Unix.out_channel_of_descr outde) result; exit 0)
      | pid -> (pid, Unix.in_channel_of_descr inde, inde, outde) in

  let close_proc _ ch fd1 fd2 = 
    let result = input_binary_int ch in 
    close_in ch; Unix.close fd2 ; result in

  let rec proc_fold pknil =
    (* If array is filled with none and pknil is finished, we're done *)
    if Array.fold_right (function None -> (&&) true | _ -> (&&) false)
         procs (finishedp pknil) then pknil else
    let new_proc pkn = if finishedp pkn then None else Some (open_proc ()) in
    let kons procn pknil =
      match procs.(procn) with
        None -> (procs.(procn) <- new_proc pknil ; pknil)
      | Some (pid, inch, infd, outfd) -> (
          match Unix.waitpid [Unix.WNOHANG] pid with
              (0, _) -> (ignore (Unix.select [] [] [] 0.05) ; pknil)
            | _ -> let pknil = pkons (close_proc pid inch infd outfd) pknil in
                (procs.(procn) <- new_proc pknil ; pknil)) in
    let pknil = Mu.fold kons pknil (Mu.range 0 (nprocs - 1)) in
    proc_fold pknil in
  proc_fold pknil 

(* Execute boolean func in batches of batch on nprocs up to max times until
   errt is reached on the fraction of time func returns true.  Returns (total,
   ntrue) of runs.  batch should be high enough that running batch number of
   funcs takes a at least a few seconds to run.  errt is the proportional size
   of the 95% CI. *)
let dogged_monte_carlop label nprocs batch max errt func =
  let finishedp errt (n, ntrue) = if n > max then true else
    let q = fl ntrue /. fl n in
    (n > 0 && ((sterr n q < errt && ntrue < n && sterr n q < errt *. q) 
                 || n >= max)) in
  let batchf () = Mu.rec_n batch (fun ct -> if func () then ct+1 else ct) 0 in
  let upf_pkons ntrue (tot, old_ntrue) = 
    Printf.printf "UPD-%s\t%d\t%d\n%!" label (tot + batch) (old_ntrue + ntrue);
    (tot + batch, old_ntrue + ntrue) in
  unix_proc_fold nprocs batchf upf_pkons (0, 0) (finishedp errt)

let pmax (naa, naA, nAA) = max (max naa naA) nAA
let ptype pop = match pop with
    (0, 0, _) -> G_AA
  | (0, _, 0) -> G_aA
  | (_, 0, 0) -> G_aa
  | _ -> failwith "ptype called when multiple types present"

let simulate nchoice wij frac pop_size init_aa init_aA =
  let stopp pop = pmax pop = pop_size || popct pop G_aa = popct pop G_AA in
  let final_pop = Mu.rec_p stopp (mch pop_size nchoice wij frac)
    (init_aa, init_aA, pop_size - init_aa - init_aA) in
  if popct final_pop G_aa = popct final_pop G_AA then Rand.float 1.0 > 0.5
  else ptype final_pop = G_aa

let (* fixation_mcmc *) () =
  let run wij f (init_aa, init_aA) nch =
    let nn = 1000 in
    let tag = Printf.sprintf "%d\t%s\t%f\t%d\t%d" 
      nn (vprf wij) f init_aa init_aA in
    let (tot, naa) = dogged_monte_carlop tag 24 10_000 100_000_000 0.1
      (fun () ->
        ignore (Rand.init (int_of_float (Unix.gettimeofday () *. 1000000.))) ;
        simulate nch wij f nn init_aa init_aA) in
    let q = fl naa /. fl tot in
    Printf.printf "%d\t%d\t%s\t%f\t%f\t%d\t%d\t%f\t%d\t%d\t%e\t%e\n%!"
      nch nn (vprf wij) (1. -. vmin wij) f init_aa init_aA
      ((fl init_aa +. 0.5 *. fl init_aA) /. fl nn)
      naa tot q (sterr tot q) in
    
  Printf.printf
    "cs\tpopsize\twij\ts\tf\tinithom\tinithet\tp\tnaa\tntrial\tq\terr\n%!" ;

  ignore (Mu.cmap4 run
    [(1.,1.,1.);
     (1.,0.999,1.);   (1.,0.99,1.);  (1.,0.9,1.); (1.,0.5,1.);
     (1.,1.,0.999);   (1.,1.,0.99);  (1.,1.,0.9); (1.,1.,0.5);
     (1.,0.999,0.999); (1.,0.99,0.99); (1.,0.9,0.9);(1.,0.5,0.5);
     (0.999,1.,1.);   (0.99,1.,1.);  (0.9,1.,1.); (0.5,1.,1.);
     (0.999,0.999,1.); (0.99,0.99,1.); (0.9,0.9,1.); (0.5,0.5,1.);
     (0.999,1.,0.999); (0.99,1.,0.99); (0.9,1.,0.9); (0.5,1.,0.5)]
    [0.;0.001;0.01;0.1;0.5;1.0]
    [(0,1);(1,2)]
    [1;2;3;5;8;13;21;34;55;89;144;233;377;610;987;1597;2584;4181;6765])
