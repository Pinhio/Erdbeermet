from typing import Union
from erdbeermet.simulation import simulate
from erdbeermet.recognition import recognize
from erdbeermet.visualize.BoxGraphVis import plot_box_graph
from time import time

def low_performing_pipeline(size:Union[int,list], iterations:int=1, first_four_simulation:list=[0,1,2,3], circular:bool=False, clocklike:bool=False, first_candidate_only:bool=True, print_info:bool=False, print_pipe_info:bool=False):
    '''
    pipeline as described in WP1

    size: simulation matrix size (list or int)
    iterations: number of iterations performed for each passed size
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

    fn = time()

    with open(f'prak/sim_outputs/{subfolder}/{fn}.txt', 'w') as f:

        for s in size:
            # measure runtime
            runtimes = []

            divergences = []

            for i in range(iterations):
                # start_time
                start_time = time()

                # create scenario 
                scenario = simulate(s, circular, clocklike)

                # get all triples of simulation
                hist = scenario.history
                simulation_triples = []
                for h in hist:
                    simulation_triples.append((h[:3]))   

                # recognition
                rec_tree = recognize(scenario.D, first_candidate_only, print_info) # additional params?
                
                # recognized as valid R-map
                rec_as_r_map = rec_tree.root.valid_ways > 0

                # reconstruction failed: print matrix, save tree, show box-graphs
                if not rec_as_r_map:
                    print("Recognition failed for matrix:")
                    f.write('===========================================\n')
                    f.write(f'recognition failed for matrix with size {s}:\n')
                    f.write(str(rec_tree.root.D))
                    f.write('\n')
                    f.write(f'history of scenario:\n')
                    f.write(str(hist))
                    f.write('\n')
                    print(rec_tree.root.D)
                    dest = f"prak/fail_outputs/{fn}"
                    rec_tree.visualize(save_as=dest + ".pdf", popup=False)
                    """ for node in rec_tree.preorder():
                        if node.valid_ways == 1 and len(node.V) == 4:
                            pass """
                            # plot_box_graph(node.D, labels=range(4))
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
                rec_first_four = first_four_simulation in valid_4_leaves

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

                # print info of this iteration
                if print_pipe_info:
                    # write to file
                    f.write(f'----------- size {s} | iteration {i+1} -----------\n')
                    # f.write(scenario.D)
                    f.write(f'simulation triples: {simulation_triples}\n')
                    f.write(f'Valid leaves : {valid_4_leaves}\n')
                    f.write(f'Valid triples: {valid_triples}\n')
                    if rec_first_four:
                        f.write(f'first four leaves ({first_four_simulation}) recognized\n')
                    f.write(f'Number of common triples: {common_triples}\n')
                    f.write(f'divergence: {divergence}\n')
                    f.write('\n')

            
            # print avg runtime of this size
            avg_runtime = sum(runtimes) / len(runtimes)
            avg_div = sum(divergences) / len(divergences)
            if print_pipe_info:
                print(f'avg runtime size {s}: {avg_runtime}')
            f.write('\n')
            f.write(f'=== size {s} ================================\n')
            f.write(f'avg runtime:    {avg_runtime}\n')
            f.write(f'avg divergence: {avg_div}\n')
            f.write('===========================================\n')
            f.write('\n')
                
            
low_performing_pipeline([6,7,8,9,10], iterations=30, print_pipe_info=False)
# other examples
# low_performing_pipeline(7, iterations=2, clocklike=True)
# low_performing_pipeline([6, 8, 10], iterations=3)