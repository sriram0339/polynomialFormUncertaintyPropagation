//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#include <deque>
#include "ProbabilityQueryEvaluator.hh"
#include <cmath>

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
                    std::cout << "}" << std::endl;
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
        polynomialExpectation = mp.expectation(distributionInfo);
        MpfiWrapper sumOfExpectations(0.0);
        for (auto & c: splitComponents){
            MpfiWrapper cExpect = c.expectation(distributionInfo);
            componentExpectations.push_back(cExpect);
            componentRanges.push_back(c.evaluate(distributionRanges));
            c.addToConst(-1.0 * cExpect);
            sumOfExpectations = sumOfExpectations + cExpect;
            MultivariatePoly cSq = c.squarePoly();
            MpfiWrapper tmp = cSq.expectation(distributionInfo);
            if (tmp.lower() < 0){
                std::cout << "Warning: second moment is negative???" << tmp << std::endl;
                c.prettyPrint(std::cout, std::map<int, string>());
                cSq.prettyPrint(std::cout, std::map<int, string>());
            }
            componentVariances.push_back(tmp);
        }
        if (debug ){
            std::cout << "Separate Polynomial into components:" << std::endl;
            std::cout << "Original Expectation: " << polynomialExpectation << std::endl;
            std::cout << "Sum of component expectations: " << sumOfExpectations << std::endl;
        }

    }

    double ProbabilityQueryEvaluator::computeChebyshevBounds() const {
        /*
         * 1. Assume that the polynomial has been separated into components.
         * 2. Sum up the variances of each component.
         */
        MpfiWrapper variancesSum (0.0);
        for (auto c: componentVariances) {
            variancesSum = variancesSum + c;
        }
        double bound = variancesSum.upper() / ( (t - polynomialExpectation.upper())* (t - polynomialExpectation.upper()) );
        return bound;
    }

    double ProbabilityQueryEvaluator::computeChernoffBound() const {
        double s = t - polynomialExpectation.upper();
        double denSum = 0;
        for (auto c: splitComponents){
            MpfiWrapper rng = c.evaluate(distributionRanges);
            denSum = denSum + (rng.upper() - rng.lower()) * (rng.upper() - rng.lower());
        }
        double expFactor = - 2* s * s / denSum;
        return exp(expFactor);
    }


    double ProbabilityQueryEvaluator::computeFourthMomentBound() const {
        MpfiWrapper fourthMomentSum(0.0);
        for (auto c: splitComponents){
            MultivariatePoly c4 = c.powPoly(4);
            MpfiWrapper fourthMoment = c4.expectation(distributionInfo);
            fourthMomentSum = fourthMoment + fourthMomentSum;
        }
        int n = componentVariances.size();
        int i,j;
        for (i =0 ; i < n; ++i){
            for (j = i+1; j < n; ++j){
                fourthMomentSum = fourthMomentSum + 2 * componentVariances[i] * componentVariances[j];
            }
        }
        double denom = (t - polynomialExpectation.upper())* (t - polynomialExpectation.upper()) *  (t - polynomialExpectation.upper()) * (t - polynomialExpectation.upper()) ;
        double bound = fourthMomentSum.upper() / denom;
        return bound;

    }


};