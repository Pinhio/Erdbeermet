from typing import Union
from erdbeermet.simulation import simulate
from erdbeermet.recognition import recognize
from erdbeermet.visualize.BoxGraphVis import plot_box_graph
from time import time
from itertools import permutations


def low_performing_pipeline(size:Union[int,list], iterations:int=1, first_four_simulation:list=[0,1,2,3], circular:bool=False, clocklike:bool=False, first_candidate_only:bool=True, block_leaves:int=0, choose_smallest_spike:bool=False, pdf_error:bool=False, print_failed:bool=False, print_info:bool=False):
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
                    if circle:
                        circles += 1
                    if print_failed:
                        f.write('===========================================\n')
                        f.write(f'recognition failed ({fails}) for matrix with size {s}:\n')
                        f.write(str(rec_tree.root.D))
                        f.write('\n')
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

            
            # print avg runtime of this size
            total_runtime = sum(runtimes)
            avg_runtime = total_runtime / len(runtimes)
            avg_div = sum(divergences) / len(divergences)
            
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
                

# tests for WP1 with error output
low_performing_pipeline(8, iterations=20000, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, circular=True, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, clocklike=True, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, circular=True, clocklike=True, print_failed=True, pdf_error=True)

# test for WP3 with error output
low_performing_pipeline(8, iterations=20000, block_leaves=4, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, circular=True, block_leaves=4, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, clocklike=True, block_leaves=4, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, circular=True, clocklike=True, block_leaves=4, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, block_leaves=3, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, circular=True, block_leaves=3, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, clocklike=True, block_leaves=3, print_failed=True, pdf_error=True)
low_performing_pipeline(8, iterations=20000, circular=True, clocklike=True, block_leaves=3, print_failed=True, pdf_error=True)

# tests for WP4 print_failed=False in order to omit buffer problems
low_performing_pipeline(8, iterations=20000, choose_smallest_spike=True)
low_performing_pipeline(8, iterations=20000, choose_smallest_spike=True, circular=True)
low_performing_pipeline(8, iterations=20000, choose_smallest_spike=True, clocklike=True)
low_performing_pipeline(8, iterations=20000, choose_smallest_spike=True, circular=True, clocklike=True)
# WP4 with error output (without generating the pdfs)
low_performing_pipeline(8, iterations=100, choose_smallest_spike=True, print_failed=True)
low_performing_pipeline(8, iterations=100, choose_smallest_spike=True, circular=True, print_failed=True)
low_performing_pipeline(8, iterations=100, choose_smallest_spike=True, clocklike=True, print_failed=True)
low_performing_pipeline(8, iterations=100, choose_smallest_spike=True, circular=True, clocklike=True, print_failed=True)

# low_performing_pipeline(8, iterations=100, choose_smallest_spike=False, pdf_error=False, print_failed=False)