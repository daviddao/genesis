/**
 * @brief
 *
 * @file
 * @ingroup test
 */

#include "plausibility/plausibility.hpp"
#include "plausibility/RMQ_succinct.hpp"

#include "tree/bipartitions.hpp"
#include "tree/newick_processor.hpp"
#include "utils/logging.hpp"

namespace genesis {

// =============================================================================
//     Plausibility
// =============================================================================

void Plausibility::SpiderpigFunction (
    const std::string& reference_tree_file,
    const std::string& small_tree_file
) {
    // read trees from files
    PlausibilityTree small_tree;
    NewickProcessor::FromFile(reference_tree_file, reference_tree);
    NewickProcessor::FromFile(small_tree_file, small_tree);

    // create preorder ids for every node, and a reference map for leaves
    std::map<std::string, int> reference_map;
    int c = 0;
    for (
        auto it = reference_tree.BeginPreorder();
        it != reference_tree.EndPreorder();
        ++it
    ) {
        it.Node()->preorder_id = c;
        if (it.Node()->IsLeaf()) {
            reference_map[it.Node()->name] = c;
        }
        ++c;
    }

    // do euler traversal and collect ids on the fly
    std::vector<int> euler_ids;
    euler_ids.reserve(4 * reference_tree.NodeCount());

    // create a map to see which index a specific preorder_id appears for the first time
    std::map<int, int> preorder_to_euler_map;

    c = 0;
    for (
        auto it = reference_tree.BeginEulertour();
        it != reference_tree.EndEulertour();
        ++it
    ) { 
        int p_id = it.Node()->preorder_id;
        euler_ids.push_back(p_id);

        if(preorder_to_euler_map.count(p_id) == 1) {
            // preorder id already exist, do nothing
        } else {
            // save the euler index
            preorder_to_euler_map[p_id] = c;
        };
        //ausgabe += std::to_string(it.Node()->preorder_id) + " ";
        ++c;
    }

    // TODO fertig machen!

    // simple trick to convert vector into array and create the RMQ datastructure
    int* euler_array = &euler_ids[0];
    auto rmq = RMQ_succinct(euler_array, 4 * reference_tree.NodeCount() - 5);
    
    // RMQ query
    //LOG_DBG << "euler[" << rmq.query(1,10000) << "]" << " is " << euler_array[rmq.query(1,1000)];
    

    // TODO this should be done for each small tree!

    std::vector<pair<int,bool>> preorder_ids;
    preorder_ids.reserve(2 * small_tree.NodeCount());

    for (
        auto it = small_tree.BeginPreorder();
        it != small_tree.EndPreorder();
        ++it
    ) {
        if(it.Node()->IsLeaf()) {
            preorder_ids.push_back(std::make_pair(reference_map[it.Node()->name],1));
        }
    }

    // sort the preorder vector the first time
    std::sort(preorder_ids.begin(),preorder_ids.end());

    // declare ints to store euler ids once for the loop
    int euler_idx,
        euler_idy,
        euler_res;

    for (
        size_t i = 0;
        i < small_tree.LeafCount() - 1;
        ++i 
    ) { 
        // get the euler index of the preorderids
        euler_idx = preorder_to_euler_map[preorder_ids[i].first];
        euler_idy = preorder_to_euler_map[preorder_ids[i + 1].first];

        // do a rmq query and get the euler id of the minimum preorder id
        euler_res = rmq.query(euler_idx,euler_idy);
        // push the preorder id into our vector
        preorder_ids.push_back(std::make_pair(euler_array[euler_res],0));
    }

    // assertion that the loop went correct
    assert(preorder_ids.size() == (2*small_tree.LeafCount() - 1));

    // sort the preorder vector the second time
    std::sort(preorder_ids.begin(),preorder_ids.end());

    // TODO: get the bipartitions from the sorted array 
    for (
        size_t i = 0;
        i < preorder_ids.size();
        ++i
    ) {
        //LOG_DBG << preorder_ids[i].first << ':' << preorder_ids[i].second << '\n';
    }


}

} // namespace genesis
