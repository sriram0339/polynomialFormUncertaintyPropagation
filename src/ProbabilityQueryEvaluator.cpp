//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#include <deque>
#include "ProbabilityQueryEvaluator.hh"

namespace PolynomialForms{
    extern bool debug;

    class NoiseSymbolGraph {
    protected:
        int numVerts;
        std::vector< std::set<int> > adjList;
        MultivariatePoly const & mp;
        std::vector< std::set<int> > components;
        std::vector<MultivariatePoly> & splitPolys;
        /*- Add an edge between i and j -*/
        void addEdge(int i, int j) {
            assert( i >= 0 && i < numVerts);
            assert( j >= 0 && j < numVerts);
            adjList[i].insert(j);
            adjList[j].insert(i);
        }

        /*- Add edges corresponding to symbols that co-occur in the same powerproduct -*/
        void addEdges(){
            std::set<PowerProduct> const & pp = mp.getPowerProducts();
            for (PowerProduct const& p: pp){
                set<int> const & indices = p.getIndices();
                for (int i: indices){
                    for (int j: indices) {
                        if (i < j){
                            addEdge(i, j);
                        }
                    }
                }
            }
        }

        /* - Do a BFS of the graph starting from a starting vertex.
         *   The method adds to the set of visited nodes and the set of nodes
         *   in the component belonging to the startVertex.
         -*/
        void doBFS(int startVert, set<int> & visitedNodes, set<int> & curComponent){
            // Insert start vertex into visited node and current component
            visitedNodes.insert(startVert);
            curComponent.insert(startVert);
            std::deque<int> queue;
            // Insert start vertex into the queue
            queue.push_front(startVert);
            while (queue.size() > 0){
                // Take the first element of the queue
                int v = queue.front();
                queue.pop_front();
                assert(v >= 0 && v < numVerts);
                // Iterate through its adjacent vertices
                set<int> const & adjVerts= adjList[v];
                for (auto j: adjVerts){
                    // If it has not been seen before
                    if (visitedNodes.find(j) == visitedNodes.end()){
                        curComponent.insert(j); // add it to the current component
                        visitedNodes.insert(j); // add it to the set of visited nodes
                        queue.push_back(j); // push it back into the queue
                    }
                }
            }
        }

        /*- Identify SCCs using DFS.
         *  This method populates the list of components.
         -*/
        void doBFSForIdentifyingSCCs(){
            set<int> visitedNodes;
            for (int i = 0; i < numVerts; ++i){
                if (visitedNodes.find(i) == visitedNodes.end()){
                    set<int> curComponent;
                    doBFS(i, visitedNodes, curComponent);
                    components.push_back(curComponent);
                }
            }
            if (debug){
                std::cout << "DEBUG: components of the polynomial are " << std::endl;
                for (auto c: components){
                    std::cout << "\t {";
                    for (int j: c){
                        std::cout << "w" << j << " ";
                    }
                    std::cout << "}";
                }
            }
        }

        /* Split the polynomial according to the components */
        void splitPolynomialAccordingToComponents(){
            if (debug){
                std::cout << "Original poly:";
                mp.prettyPrint(std::cout, std::map<int, string>());
                std::cout << std::endl;
            }
            // Iterate through the components
            std::map<PowerProduct, MpfiWrapper> const & terms = mp.getTerms();

            for (set<int> const & c: components){
                MultivariatePoly p(0.0); // New polynomial
                for (auto t: terms){
                    // Iterate through power products
                    PowerProduct const & pp = t.first;
                    if (pp.hasIntersectionWith(c)){
                        p.setTerm(t.first, t.second);
                    }
                }
                splitPolys.push_back(p);
                if (debug){
                    std::cout << "Split poly: " ;
                    p.prettyPrint(std::cout, std::map<int, string>());
                    std::cout << std::endl;
                }
            }
        }

    public:
        NoiseSymbolGraph(int n, MultivariatePoly const & mp_, std::vector<MultivariatePoly> & res_): numVerts(n),
        adjList(n, std::set<int>()), mp(mp_), splitPolys(res_) {
            addEdges();
            doBFSForIdentifyingSCCs();
            splitPolynomialAccordingToComponents();
        }

    };


    void ProbabilityQueryEvaluator::separatePolynomialIntoComponents() {
        NoiseSymbolGraph nsGraph(distributionInfo.size(), mp, splitComponents);
    }


};