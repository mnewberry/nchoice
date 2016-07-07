module M = Moran
module Pop = M.Pop
module Meas = M.Meas
let fl = float_of_int

type genotype = G_aa | G_aA | G_AA

let gen_prob dist offsp = match offsp, dist with
      G_aa, (p,_,_) -> p
    | G_aA, (_,p,_) -> p
    | G_AA, (_,_,p) -> p

let mendel_dist a b = match (a,b) with
      G_aa, G_aa -> (1., 0., 0.)
    | G_aa, G_aA -> (0.5, 0.5, 0.)
    | G_aa, G_AA -> (0., 1., 0.)
    | G_aA, G_aa -> (0.5, 0.5, 0.)
    | G_aA, G_aA -> (0.25, 0.5, 0.25)
    | G_aA, G_AA -> (0., 0.5, 0.5)
    | G_AA, G_aa -> (0., 1., 0.)
    | G_AA, G_aA -> (0., 0.5, 0.5)
    | G_AA, G_AA -> (0., 0., 1.)

let mendel_prob offsp a b = gen_prob (mendel_dist a b) offsp

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

let update selfp s ((p_aa, p_aA, p_AA) as initps) =
  let lu = gen_prob initps in
  let genotypes = [G_aa;G_aA;G_AA] in
  let prob zygote =
    let pop_d = (lu G_aa, lu G_aA, lu G_AA) in
    let mating_prob i j = let p = gen_prob pop_d in
      if i = j then selfp *. p i +. (1. -. selfp) *. p i *. p i
      else (1. -. selfp) *. p i *. p j in
    Mu.sumf (Mu.cmap2
              (fun i j -> mating_prob i j *. mendel_prob zygote i j)
              genotypes genotypes) in
  let pp_aa, pp_aA, pp_AA = (prob G_aa, prob G_aA, prob G_AA) in
  let ec_aa, ec_aA, ec_AA = (pp_aa, (1. -. s) *. pp_aA, pp_AA) in
  let norm x = x /. (ec_aa +. ec_aA +. ec_AA) in
  (norm ec_aa, norm ec_aA, norm ec_AA)

let freqv_of_pop p = 
  let lu g = fl (try Pop.lookup g p with Not_found -> 0)/.fl (Pop.size p) in
  (lu G_aa, lu G_aA, lu G_AA)

let meas_of_freqv (paa, paA, pAA) =
  Mu.fold2 Meas.fadd Meas.empty [paa; paA; pAA] [G_aa; G_aA; G_AA]

let pop_of_ps_aa_aA ps naa naA =
  Mu.fold2 Pop.add Pop.empty [(ps - naa - naA); naA; naa] [G_AA; G_aA; G_aa]

let simulate selfp sel pop_size init_aa init_aA =
  let wf_sample_next pop =
    let freqv = freqv_of_pop pop in
    let dist = meas_of_freqv (update selfp sel freqv) in
    Meas.sample_pop (Pop.size pop) dist in
  let stopp pop = (Pop.max_pop pop = pop_size && Pop.max_type pop <> G_aA) in
  let final_pop = Mu.rec_p stopp
    (wf_sample_next)
    (pop_of_ps_aa_aA pop_size init_aa init_aA) in
  (Pop.max_type final_pop = G_aa)

let _ =
  let run selfp nn sel () =
    let init_aa, init_aA = 0, 1 in
    let tag = Printf.sprintf "%d\t%f\t%d\t%d" nn sel init_aa init_aA in
    let (tot, naa) = dogged_monte_carlop tag 22 100_000 100_000_000 0.1
      (fun () ->
        ignore (Rand.init (int_of_float (Unix.gettimeofday () *. 1000000.))) ;
        simulate selfp sel nn init_aa init_aA) in
    let q = fl naa /. fl tot in
    Printf.printf "assortmc\t%f\t%d\t%f\t%d\t%d\t%f\t%d\t%d\t%e\t%e\n%!"
      selfp nn sel init_aa init_aA
      ((fl init_aa +. 0.5 *. fl init_aA) /. fl nn)
      naa tot q (sterr tot q) in
    
  Printf.printf
    "model\tselfp\tpopsize\ts\tf\tinithom\tinithet\tp\tnaa\tntrial\tq\terr\n%!" ;

  Mu.cmap4 run
    [0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1.0]
    [1000]
    [0.;0.001;0.01;0.1]
    [()]
