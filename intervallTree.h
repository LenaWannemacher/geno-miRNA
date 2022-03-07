#ifndef MIRNA_INTERVALLTREE_H
#define MIRNA_INTERVALLTREE_H

#include <iostream>

struct Interval{

    int low{}, high{};
    std::string loc;
    //std::vector<Interval> gene_region;

};

struct ITNode{

    Interval *i;
    int max;
    ITNode *left, *right;
    ITNode *gene_root;

};

/// 25.10.

struct Gene{

    // A gene is represented as an Interval (the genes borders) and a root of an IT, that contains all the genes intern borders and locations
    Interval gene_interval{};
    ITNode *gene_list_node{};

};


inline ITNode * newNode(const Interval& i, ITNode *gene_root_node){
    auto *temp = new ITNode;
    temp->i = new Interval(i);
    temp->max = i.high;
    temp->left = temp->right = NULL;
    temp->gene_root = gene_root_node;
    return temp;
}

inline ITNode * newNode(const Interval& i){

    auto *temp = new ITNode;
    temp->i = new Interval(i);
    temp->max = i.high;
    temp->left = temp->right = NULL;
    temp->gene_root = NULL;
    return temp;

}

/**
 * inserts a new node into the IT
 *
 * @param root  The given root_node
 * @param i     The interval corresponding to the node to be inserted
 * @param gene_root     The corresponding gene_root_node
 * @return  Pointer to the new node
 */
inline ITNode *insert_with_node(ITNode *root, Interval i, ITNode *gene_root){

    if(root == NULL){
        return newNode(i, gene_root);
    }

    int l = root->i->low;

    if(i.low < l){
        root->left = insert_with_node(root->left, i, gene_root);
    }
    else {
        root->right = insert_with_node(root->right, i, gene_root);
    }

    if(root->max < i.high){
        root->max = i.high;
    }

    return root;
}

/**
 * inserts a new node into the IT
 *
 * @param root  The given root_node
 * @param i     The interval corresponding to the node to be inserted
 * @return  Pointer to the new node
 */
inline ITNode *insert(ITNode *root, Interval i){

    if(root == NULL){
        return newNode(i);
    }

    int l = root->i->low;

    if(i.low < l){
        root->left = insert(root->left,i);
    }
    else{
        root->right = insert(root->right, i);
    }

    if(root->max < i.high){
        root->max = i.high;
    }

    return root;

}

/**
 * checks if two intervals overlap
 *
 * @param i1    The first interval
 * @param i2    The second interval
 * @return  True if they overlap, otherwise false
 */
inline bool doOverlap(Interval i1, Interval i2){

    if(i1.low <= i2.high && i2.low <= i1.high){
        return true;
    }
    return false;
}

/**
 * recursive function that searches for the node in an IT that achieves the best overlap with a given interval
 *
 * @param root  The root_node of the IT of interest
 * @param i     The interval that we search the overlap for
 * @return  The interval that achieved the best overlap
 */
inline Interval *overlapSearch(ITNode *root, Interval i){

    if(root == NULL) {
        return NULL;
    }

    if(doOverlap(*(root->i),i)){
        return root->i;
    }

    if(root->left != NULL && root->left->max >= i.low){
        return overlapSearch(root->left, i);
    }

    return overlapSearch(root->right, i);
}
/**
 * Function that processes an IT for the overlap search
 *
 * @param root  The root_node of the IT of interest
 * @param i     The interval that we search the overlap for
 * @param poss_hits     An interval vector that saves all overlaps. Can be used to determine the best possible overlap
 */
inline void process_for_overlapSearch(ITNode *root, Interval i, std::vector<Interval>& poss_hits){

    if(root == nullptr){
        return;
    }

    if(doOverlap(*(root)->i, i)){
        poss_hits.push_back(*(root)->i);
    }

    if(root->left != nullptr && root->left->max >= i.low){
        process_for_overlapSearch(root->left, i, poss_hits);
    }

    if(root->right != nullptr) {
        process_for_overlapSearch(root->right, i, poss_hits);
    }

}

/**
 * Function that processes an IT for the overlap search
 *
 * @param root  The root_node of the IT of interest
 * @param i     The interval that we search the overlap for
 * @param poss_node_hits    An ITNode* vector that saves all overlaps. Can be used to determine the best possible overlap
 */
inline void process_for_overlapSearch_with_nodes(ITNode *root, Interval i, std::vector<ITNode*>& poss_node_hits){

    if(root == nullptr){ return; }

    if(doOverlap(*(root)->i, i)){
        poss_node_hits.push_back(root);
    }

    if(root->left != nullptr && root->left->max >= i.low){
        process_for_overlapSearch_with_nodes(root->left, i, poss_node_hits);
    }

    if(root->right != nullptr){ //&& root->right->i->low <= i.low) {
        process_for_overlapSearch_with_nodes(root->right, i, poss_node_hits);
    }

}

/**
 * Finds the best possible overlap
 *
 * @param poss_hits An interval vector containing all intervals that overlapped with the interval of interest
 * @param i     The interval of interest
 * @return  The interval with the best possible overlap
 */
inline Interval best_Overlap(std::vector<Interval>& poss_hits, Interval& i){

    Interval best_int;

    if(poss_hits.empty()){
        Interval no_hit = {0,0,"no_hit"};
        return no_hit;
    }

    int best, current;
    int left, right;

    bool first = true;

    for(const auto& it: poss_hits){

        left = it.low - i.low;
        right = it.high - i.high;

        if(left <= 0 && right >= 0){

            current = abs(left) + right;

            if(first){
                best_int = it;
                best = current;
                first = false;
            }
            else{
                if(current <= best){
                    best = current;
                    best_int = it;
                }
            }
        }
        else {
            best_int.loc = "no_hit";
        }
    }

    return best_int;

}

/**
 * Prints an IT in order
 *
 * @param root  The root_node of the IT
 */
inline void inorder(ITNode *root){

    if(root == NULL){ return; }

    inorder(root->left);

    std::cout << "[" << root->i->low << ", " << root->i->high << "]" << " max = " << root->max << std::endl;

    inorder(root->right);

}


#endif //MIRNA_INTERVALLTREE_H
