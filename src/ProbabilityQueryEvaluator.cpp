//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#include <deque>
#include "ProbabilityQueryEvaluator.hh"
#include <cmath>

namespace PolynomialForms{
    extern bool debug;
    extern bool fourthMomentBoundCalculation;

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
        void splitPolynomialAccordingToComponents_cutEdges(map<int, MpfiWrapper> const & distributionRanges,MpfiWrapper & interval_cut){
            std::vector< std::set<int> > components_cut;
            std::vector<MultivariatePoly>  splitPolys_cut;
            set<PowerProduct> CutTerms;
            for (int i=0;i<components.size();i++){
                set<int> c=components[i];
                int nc=c.size();
                if(nc >=4){
                    map<int,int> mp_c;//map the original index of current component into new index from 0-(nc-1)
                    map<int,int> mp_cb;//map the new index back to the original index
                    int count=0;
                    for(auto & c_temp:c){
                        mp_c[c_temp]=count;
                        mp_cb[count]=c_temp;
                        count++;
                    }
                    vector<vector<double>> w(nc,vector<double> (nc));
                    //compute the weight of edges
                    std::map<PowerProduct, MpfiWrapper> const & terms=splitPolys[i].getTerms();
                    for(auto t: terms){
                        PowerProduct const & pp = t.first;
                        set<int> const & indices = pp.getIndices();
                        double absM=fabs(median(t.second));
                        for(auto & id1:indices){
                            int newId1=mp_c[id1];
                            for(auto & id2:indices){
                                int newId2=mp_c[id2];
                                if(newId1!=newId2){
                                    w[newId1][newId2]+=absM;
                                }
                            }
                        }
                    }
                    vector <double> ew;
                    vector <pair<int,int>> edges;
                    for(int i=0;i<nc;i++){
                        for(int j=0;j<i;j++){
                            if(w[i][j]!=0){
                                ew.push_back(w[i][j]);
                                edges.push_back(make_pair(i,j));
                            }
                        }
                    }
                    int K=3;//ceil(double(nc)/2.0);
                    vector<vector<int>> Clusters(K,vector<int>());
                    partitioning(ew,edges,nc,K,Clusters);
                    for(int k=0;k<K;k++){
                        vector<int> const & cluster_k = Clusters[k];
                        set<int> c_k;
                        //map to original index
                        for(auto a:cluster_k){
                            c_k.insert(mp_cb[a-1]);//index in glpk:1-nc, index in map: 0-(nc-1)
                        }
                        components_cut.push_back(c_k);
                        MultivariatePoly p(0.0); // New polynomial
                        for (auto t: terms) {
                            // Iterate through power products
                            PowerProduct const &pp = t.first;
                            set<int> const &indices = pp.getIndices();
                            /* Check if all of the vars in the power product are in the set vars of c_k */
                            bool keepTerm = std::all_of(indices.begin(), indices.end(), [&](int i) {
                                return c_k.find(i) != c_k.end();
                            });
                            if (pp.hasIntersectionWith(c_k)){
                                if(keepTerm) {
                                    p.setTerm(t.first, t.second);
                                }else if(!keepTerm && !CutTerms.count(pp)){
                                    MpfiWrapper interval_t = t.second * t.first.evaluate(distributionRanges);
                                    interval_cut = interval_cut + interval_t;
                                    CutTerms.insert(pp);
                                }
                            }
                        }
                        splitPolys_cut.push_back(p);
                        if (debug){
                            std::cout << "Split poly by cutting edges: " ;
                            p.prettyPrint(std::cout, std::map<int, string>());
                            std::cout << std::endl;
                        }
                        
                    }
                    
                }else{
                    components_cut.push_back(components[i]);
                    splitPolys_cut.push_back(splitPolys[i]);
                }
            }
            //            if(splitPolys_cut.size()>0){
            //                splitPolys_cut[0].addToConst(interval_cut);
            //            }
            splitPolys=splitPolys_cut;
            components=components_cut;
        }
        void partitioning(vector<double> const & w,vector <pair<int,int>> const & edges,int d,int K,vector<vector<int>> & Clusters)
        {
            //edges{vi,vj}, index i,j=0.,1..
            //w:weigh of eij
            //int d:number of state variables
            //int K=2;//Clusters_num
            int L=1;//at least L nodes in each cluster
            
            //vertics:d  + K
            //integer variables:v_ik, e_ij
            //[v_11 v_12 ... v_1K v_21 v_22 ... v_2K...v_d1 v_d2 ... v_dK e_21 e_31 e_32...] where eij is non-zero
            int ne=edges.size();
            int nx=d*K+ne;//number of variables
            int nrow=d+K+K*ne;//number of constraints
            int size = 2*d*K+3*K*ne;
            
            int *rowInd = new int[ 1 + size ];
            int *colInd = new int[ 1 + size ];
            double *coes = new double [ 1 + size ];
            
            glp_prob *mip = glp_create_prob();
            glp_set_obj_dir(mip, GLP_MIN);
            glp_term_out(GLP_OFF);
            
            glp_add_rows(mip, nrow);// define how many constraints
            for(int i=1; i<=d; ++i)
            {   //define the coefficient of left-hand side in (bound of) constraint; here GLP_FX:equal bound.
                glp_set_row_bnds(mip, i, GLP_FX, 1, 1);
            }
            for(int i=d+1; i<=d+K; ++i){
                glp_set_row_bnds(mip, i, GLP_LO, L, L);
            }
            for(int i=d+K+1; i<=nrow; ++i){
                glp_set_row_bnds(mip, i, GLP_UP, 0, 0);
            }
            
            glp_add_cols(mip, nx);// define how many variables
            for(int i=1; i<=nx; ++i){
                glp_set_col_bnds(mip, i, GLP_DB, 0, 1);
                glp_set_col_kind(mip, i, GLP_IV);
            }
            //set the objective: %e_ij
            for(int i=d*K+1; i<=nx; ++i){
                glp_set_obj_coef(mip, i, w[i-d*K-1]);
            }
            
            //set the matrix of constraints
            int pos=1;
            for(int i=1; i<=d; i++){
                for(int j=(i-1)*K+1; j<=i*K; j++){
                    //int pos = j + (i-1)*nx;
                    rowInd[pos] = i;
                    colInd[pos] = j;
                    coes[pos] = 1;
                    pos++;
                }
            }
            for(int i=1; i<=K; i++){
                for(int j1=1; j1<=d; j1++){
                    int j=(j1-1)*K+i;
                    //int pos = j + (i-1)*nx +d*nx;
                    rowInd[pos] = i+d;
                    colInd[pos] = j;
                    coes[pos] = 1;
                    pos++;
                }
            }
            for(int i1=1; i1<=ne; i1++){
                for(int k=1; k<=K; k++) {
                    int index_j=d*K+i1;//e_ij
                    rowInd[pos] = (i1-1)*K+k + (d+K);
                    colInd[pos] = index_j;
                    coes[pos] = -1;
                    pos++;
                    int i = edges[i1-1].first + 1;
                    int j = edges[i1-1].second + 1;
                    index_j=(i-1)*K+k;//v_ik
                    rowInd[pos] = (i1-1)*K+k + (d+K);
                    colInd[pos] = index_j;
                    coes[pos] = 1;
                    pos++;
                    index_j=(j-1)*K+k;//v_jk
                    rowInd[pos] = (i1-1)*K+k + (d+K);
                    colInd[pos] = index_j;
                    coes[pos] = -1;
                    pos++;
                }
            }
            
            //            for(int i=1;i<=size;i++){
            //                cout<<rowInd[i]<<", "<<colInd[i]<<", "<<coes[i]<<endl;
            //            }
            
            glp_load_matrix(mip, size, rowInd, colInd, coes);
            
            glp_iocp parm;
            glp_init_iocp(&parm);
            parm.presolve = GLP_ON;
            int err = glp_intopt(mip, &parm);
            if(err != 0){
                cout<< "solving MIP problem: Failed"<<endl;
                return;
            }
            
            //            double z = glp_mip_obj_val(mip);
            ////            for(int i=1;i<=nx;i++){
            ////                cout<<glp_mip_col_val(mip, i)<<endl;
            ////            }
            //            cout<<z<<endl;
            
            for(int i=1;i<=d;i++){
                for(int k=1;k<=K;k++) {
                    if(glp_mip_col_val(mip, (i-1)*K+k)==1){
                        Clusters[k-1].push_back(i);
                    }
                }
            }
            //            for(int k=0;k<K;k++) {
            //                for(int j=0;j<Clusters[k].size();j++) {
            //                    cout<<Clusters[k][j]<<" ";
            //                }
            //                cout<<endl;
            //            }
            
            glp_delete_prob(mip);
            delete[] rowInd;
            delete[] colInd;
            delete[] coes;
            
        }

    public:
        NoiseSymbolGraph(int n, MultivariatePoly const & mp_, std::vector<MultivariatePoly> & res_): numVerts(n),
        adjList(n, std::set<int>()), mp(mp_), splitPolys(res_) {
            addEdges();
            doBFSForIdentifyingSCCs();
            splitPolynomialAccordingToComponents();
        }
        NoiseSymbolGraph(map<int, MpfiWrapper> const & distributionRanges_, int n, MultivariatePoly const & mp_, std::vector<MultivariatePoly> & res_, MpfiWrapper & interval_cut_): numVerts(n),
        adjList(n, std::set<int>()), mp(mp_), splitPolys(res_){
            addEdges();
            doBFSForIdentifyingSCCs();
            splitPolynomialAccordingToComponents();
            splitPolynomialAccordingToComponents_cutEdges(distributionRanges_,interval_cut_);
        }

    };


    void ProbabilityQueryEvaluator::separatePolynomialIntoComponents() {
        //NoiseSymbolGraph nsGraph(distributionInfo.size(), mp, splitComponents);
        MpfiWrapper interval_cut(0.0);
        NoiseSymbolGraph nsGraph(distributionRanges, distributionInfo.size(), mp, splitComponents, interval_cut);
        std::cout << "Debug: Cost of interval cut = " << interval_cut << std::endl;
        t = t- interval_cut.upper(); // if p + [l,u] >= 0  then p >= -u Pr(p + [l,u] >= 0) <= P(p >= -u) <= ..
        //std::cout << "Debug: Converting p +  " << interval_cut << " >= s+t into p >= s + " << t << std::endl;

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
        if (t <= polynomialExpectation.upper()) { return 1.0;}
        MpfiWrapper variancesSum (0.0);
        for (auto c: componentVariances) {
            variancesSum = variancesSum + c;
        }
       // std::cout << variancesSum << std::endl;
        double bound = variancesSum.upper() / ( (t - polynomialExpectation.upper())* (t - polynomialExpectation.upper()) );
        return bound;
    }

    double ProbabilityQueryEvaluator::computeChernoffBound() const {
        std::cout << "Chernoff: E(X) = " << polynomialExpectation << std::endl;
        double s = t - polynomialExpectation.upper();
        std::cout << "Debug: Chernoff bounds -- X - E(X) >= " << s << std::endl;
        if ( s <= 0) { return 0.0; }
        double denSum = 0;
        for (auto c: splitComponents){
            MpfiWrapper rng = c.evaluate(distributionRanges);
            denSum = denSum + (rng.upper() - rng.lower()) * (rng.upper() - rng.lower());
        }
        std::cout << "Debug: Chernoff -- den sum = " << denSum << std::endl;
        double expFactor = - 2 * s * s / denSum;
        return expFactor;
    }


    double ProbabilityQueryEvaluator::computeFourthMomentBound() const {
        MpfiWrapper fourthMomentSum(0.0);
        if (t <= polynomialExpectation.upper()) { return 1.0;}
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


    void ProbabilityQueryEvaluator::computeBestUpperTailBounds(double s) {
        if (splitComponents.size() <= 0){
            separatePolynomialIntoComponents();
        }
        std::cout << "Evaluating probability query with bound " << s << std::endl;


        t = t + s;
        printPolyStats();
        double chBoundsLog = computeChernoffBound();
        double chBounds = exp(chBoundsLog);
        std::cout << "\t Chernoff Bounds: e^{" << chBoundsLog  << "} = " << chBounds << std::endl;
        std::cout << "\t Chebyshev Bounds:  " << computeChebyshevBounds() << std::endl;
        double bernsteinBoundsLog = computeBernsteinBound();
        double bernsteinBounds = exp(bernsteinBoundsLog);
        std::cout << "\t Bernstein Bounds: e^{" << bernsteinBoundsLog << "}=" << bernsteinBounds << std::endl;

        if (fourthMomentBoundCalculation){
            double fourthMomentBound = computeFourthMomentBound();
            std::cout << "\t Fourth Moment Bounds: " << fourthMomentBound << std::endl;
        }
        t = t - s;
    }

    double ProbabilityQueryEvaluator::computeBernsteinBound() const {
        /*
         * Bernstein bounds.
         *    X1... Xn are random variables with E(Xi) = 0
         */
        double s = t - polynomialExpectation.upper();
        if (s <= 0) { return 0.0; }
        double M = 0;
        for (auto c: splitComponents){
            MpfiWrapper crng = c.evaluate(distributionRanges);
            double Mhat = max(fabs(crng.upper()), fabs(crng.lower()));
            M = max(M, Mhat);
        }
        double sumOfVariances = 0.0;
        for (auto cvar: componentVariances) {
            sumOfVariances = sumOfVariances + cvar.upper();
        }
        std::cout << "Debug: Bernstein inequality: M = " << M << std::endl;
        double expFactor = -0.5 * s*s/(sumOfVariances + M * s/3.0 );
        return expFactor;
    }

    void ProbabilityQueryEvaluator::printPolyStats() {
        std::cout << "Stats about polynomial form in query evaluation. " << std::endl;
        std::cout << "Number of noise symbols: " << distributionInfo.size() << std::endl;
        std::cout << "Total Degree: " << mp.degree() << std::endl;
        std::cout << "Num Terms: " << mp.getTerms().size()<< std::endl;
        std::cout << "Number of split polynomials: " << splitComponents.size() <<std::endl;
        int largestSplitComp = 0;
        for (auto c: splitComponents){
            int j = c.getTerms().size();
            if (j > largestSplitComp){
                largestSplitComp = j;
            }
        }
        std::cout << "Largest size split component (number of terms): " << largestSplitComp << std::endl;
    }


};
