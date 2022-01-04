from typing import Union
from erdbeermet.simulation import simulate
from erdbeermet.recognition import recognize
from erdbeermet.visualize.BoxGraphVis import plot_box_graph
from time import time

def low_performing_pipeline(size:Union[int,list], first_four:list=None, circular:bool=False, clocklike:bool=False, iterations:int=1, first_candidate_only:bool=False):
    '''
    pipeline as described in WP1

    size: simulation matrix size (list or int)
    '''
    # create list of single size to guarantee functionality
    if type(size) == int:
        size = [size]

    # measure runtime
    runtimes = []

    for s in size:
        for _ in range(iterations):
            # start_time
            start_time = time()

            # create scenario 
            scenario = simulate(s, circular, clocklike)

            # get all triples
            hist = scenario.history
            triple = []
            for h in hist:
                triple.append((h[:3]))
            print(f'triples: {triple}')
            # print(hist)
            
            # get first four leaves if not passed as param
            if first_four == None:
                first_four = [triple[0][0]]
                for t in range(3):
                    first_four.append((triple[t][2]))
                print(f"First four leaves: {first_four}")
                

            # recognition
            rec_tree = recognize(scenario.D, first_candidate_only) # additional params?
            
            # recognized as valid R-map
            rec_as_r_map = rec_tree.root.valid_ways > 0

            # reconstruction failed: print matrix, save tree, show box-graphs
            if not rec_as_r_map:
                print("Recognition failed for matrix:")
                print(rec_tree.root.D)
                dest = "prak/fail_outputs/" + str(time())
                rec_tree.visualize(save_as=dest + ".pdf", popup=False)
                for node in rec_tree.preorder():
                    if node.valid_ways == 1 and len(node.V) == 4:
                        pass
                        # plot_box_graph(node.D, labels=range(4))
                # continue

            # get all valid combinations of 4 leaves
            valid_4_leaves = []
            valid_triples = []
            for node in rec_tree.preorder():
                if node.valid_ways == 1 and len(node.V) == 4:
                    valid_4_leaves.append(node.V)
                    valid_triples.append((node.R_step[:3]))
            print(f"Valid leaves : {valid_4_leaves}")
            print(f"Valid triples: {valid_triples}")
            

            #Classify whether the final 4-leaf map matches the first 4 leaves of the simulation
            rec_first_four = first_four in valid_4_leaves
            if rec_first_four:
                print('first four leaves recognized')

            # Count common triples
            common_triples = 0
            for t in valid_triples:
                if t in triple:
                    common_triples += 1
            divergence = common_triples / len(valid_triples)
            print(f"Number of common triples: {common_triples}")
            print(f'divergence: {divergence}')

            # end_time
            end_time = time()
            # calc runtime
            runtimes.append(end_time - start_time)


    # print(runtimes)
    avg_runtime = sum(runtimes) / len(runtimes)
    print(f"Avg runtime: {avg_runtime}")
            
            
low_performing_pipeline(6)
# other examples
# low_performing_pipeline(7, iterations=2, clocklike=True)
# low_performing_pipeline([6, 8, 10], iterations=3)