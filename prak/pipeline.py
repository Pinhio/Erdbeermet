from typing import Union
from erdbeermet.simulation import simulate
from erdbeermet.recognition import recognize
from time import time
from itertools import permutations
import os


def pipeline(size:Union[int,list], iterations:int=1, first_four_simulation:list=[0,1,2,3], circular:bool=False, clocklike:bool=False, first_candidate_only:bool=True, block_leaves:int=0, choose_smallest_spike:bool=False, generate_hist_files:bool = False, pdf_error:bool=False, print_failed:bool=False, print_info:bool=False):
    '''
    pipeline as described in WP1

    size: simulation matrix size (list or int)
    iterations: number of iterations performed for each passed size
    first_four_simulation: first four leafs
    circular: parameter passed into simulation to genereate circular matrix
    clocklike: parameter passed into simulation to genereate clocklike distances
    first_candidate_only: terminates regocnition when first valid candidate has been found
    block_leaves: initiates blocking of first 3 resp. 4 leaves of generation. These leaves won't be considered in recognition as candidates.
    choose_smallest_spike: only considers candidates with smallest spikes for recognition
    pdf_error: show pdf trees if true
    print_fail: write additional debug-information for each iteration to file
    print_info: parameter passed into recognition algorithm to print out info
    '''
    
    # create list of single size to guarantee functionality
    # don't start if size < 6 is given
    if type(size) == int:
        if size < 6:
            print('min size: 6')
            return
        size = [size]

    # find subfolder
    subfolder = '' 
    if not circular and not clocklike: subfolder = 'ncnc'
    elif circular and clocklike: subfolder = 'cc'
    elif circular and not clocklike: subfolder = 'cnc'
    elif not circular and clocklike: subfolder = 'ncc'

    # determine WP for filename
    wp = 'wp1'
    if choose_smallest_spike:
        wp = 'wp4'
    elif block_leaves in [3,4]:
        wp = 'wp3'

    # filename var which corresponds to starting time
    fn = time()
    if generate_hist_files:
        os.makedirs(os.path.dirname(f'prak/sim_outputs/{subfolder}/{wp}_{fn}_hists/'), exist_ok=True)

    with open(f'prak/sim_outputs/{subfolder}/{wp}_{fn}.txt', 'w') as f:
        # write chosen parameters to file
        f.write(f'=====================================================\n')
        f.write(f'size:                  {size}\n')
        f.write(f'iterations:            {iterations}\n')
        f.write(f'first_four_simulation: {first_four_simulation}\n')
        f.write(f'circular:              {circular}\n')
        f.write(f'clocklike:             {clocklike}\n')
        f.write(f'first_candidate_only:  {first_candidate_only}\n')
        f.write(f'block_leaves:          {block_leaves}\n')
        f.write(f'choose_smallest_spike: {choose_smallest_spike}\n')
        f.write(f'pdf_error:             {pdf_error}\n')
        f.write(f'print_failed:          {print_failed}\n')
        f.write(f'print_info:            {print_info}\n')
        f.write(f'=====================================================\n')

        for s in size:
            # measure runtime
            runtimes = []
            # measure divergence
            divergences = []
            # counter for fails
            fails = 0
            # counter for circles
            circles = 0
            # counter for recognition of the initial leaves of the simulation
            rec_first_four = 0

            for i in range(iterations):
                # start_time
                start_time = time()

                # create scenario 
                scenario = simulate(s, circular=circular, clocklike=clocklike)

                # get all triples of simulation
                hist = scenario.history
                simulation_triples = []
                for h in hist:
                    simulation_triples.append((h[:3]))   
                    
                # recognition
                B = None
                # WP3 (this if block)
                if block_leaves in [3,4]:
                    perms = permutations(range(s), block_leaves)
                    for B in perms:
                        rec_tree, circle = recognize(scenario.D, B, choose_smallest_spike, first_candidate_only, print_info)
                        if rec_tree.root.valid_ways > 0:
                            # print("valid permutation found")
                            break
                # WP4
                elif choose_smallest_spike:
                    rec_tree, circle = recognize(scenario.D, B, choose_smallest_spike, first_candidate_only, print_info)
                # WP1 (normal pipeline)
                else:
                    rec_tree, circle = recognize(scenario.D, B, choose_smallest_spike, first_candidate_only, print_info)

                if circle and pdf_error:
                    rec_tree.visualize(save_as=dest + ".pdf", popup=False)

                # recognized as valid R-map
                rec_as_r_map = rec_tree.root.valid_ways > 0

                # reconstruction failed: print matrix, save tree, show box-graphs
                if not rec_as_r_map:
                    fails += 1
                    # save hist file 
                    if generate_hist_files:
                        scenario.write_history(os.path.join(f'prak/sim_outputs/{subfolder}/{wp}_{fn}_hists/', f'error_{fails}'))   

                    if circle:
                        circles += 1
                    if print_failed:
                        f.write('===========================================\n')
                        f.write(f'recognition failed ({fails}) for matrix with size {s}:\n')
                        f.write(f'history of scenario:\n')
                        f.write(str(hist))
                        f.write('\n')
                        dest = f"prak/sim_outputs/{subfolder}/{wp}_{fn}_s{s}_{fails}"
                        if pdf_error:
                            rec_tree.visualize(save_as=dest + ".pdf", popup=False)
                    continue

                # get the first valid combination of 4 leaves
                valid_4_leaves = []
                valid_triples = []
                for node in rec_tree.preorder():
                    if node.valid_ways == 1 and len(node.V) == 4: 
                        valid_4_leaves.append(node.V)
                        valid_triples.append(node.R_step[:3])
                        break # valid result found
                    elif node.valid_ways != 0 and 4 < len(node.V) < s:
                        valid_triples.append(node.R_step[:3])

                # classify whether the final 4-leaf map matches the first 4 leaves of the simulation
                rec_first_four += 1 if first_four_simulation in valid_4_leaves else 0

                # Count common triples
                common_triples = 0
                for t in valid_triples:
                    reverse_t = (t[1], t[0], t[2]) # swap parents
                    if t in simulation_triples or reverse_t in simulation_triples:
                        common_triples += 1
                divergence = 1 - (common_triples / len(valid_triples))
                divergences.append(divergence)

                # end_time
                end_time = time()
                # calc runtime
                runtimes.append(end_time - start_time)


            # calculate avg runtime and divergence
            total_runtime = sum(runtimes)
            avg_runtime = total_runtime / len(runtimes)
            avg_div = sum(divergences) / len(divergences)

            # write kpis to file
            f.write('\n')
            f.write(f'=== size {s} ==========================================\n')
            f.write(f'total runtime:       {total_runtime}\n')
            f.write(f'avg runtime:         {avg_runtime}\n')
            f.write(f'avg divergence:      {avg_div}\n')
            f.write(f'first_four_rec:      {rec_first_four} of {iterations}\n')
            f.write(f'failed recognitions: {fails} of {iterations}\n')
            f.write(f'circles:             {circles}\n')
            f.write('=====================================================\n')
            f.write('\n')
                
# SET PARAMETERS HERE  
# int or list - max number of elements to be simulated (matrix-size)
SIZE = 8 

# int - number of simulations to carry out
ITERATIONS = 20000 

# list - specify the initial elements of the simulation
FIRST_FOUR_SIMULATION = [0,1,2,3] 

# bool - parameter for input matrix generation
CIRCULAR = True 

# bool - parameter for input matrix generation
CLOCKLIKE = True 

# bool - always choose first valid candidate
FIRST_CANDIDATE_ONLY = True 

# int - set to 3 or 4 - inhibits algorithm from choosing initial X leaves for reduction
BLOCK_LEAVES = 0

# bool - only choose smallest spike for reduction (WP4)
CHOOSE_SMALLEST_SPIKE = False 

# bool - generate and save the histories to file
GENERATE_HIST_FILES = False 

# bool - generate PDF for every failed recocnition (HANDLE WITH CARE!!!!eins!!elf!!!!)
PDF_ERROR = False 

# bool - save detailed info for failed recocnitions
PRINT_FAILED = False 

# bool - parameter for console output
PRINT_INFO = False 


# RUN PIPELINE HERE
pipeline(SIZE, ITERATIONS, FIRST_FOUR_SIMULATION, 
    CIRCULAR, CLOCKLIKE, FIRST_CANDIDATE_ONLY, BLOCK_LEAVES, 
    CHOOSE_SMALLEST_SPIKE, GENERATE_HIST_FILES, PDF_ERROR, PRINT_FAILED, PRINT_INFO)