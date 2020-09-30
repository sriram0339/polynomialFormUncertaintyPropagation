#include "inference.h"

void mono_next( int m, vector <int> & x )
{
    int n=x.size();
    x[n-1]=x[n-1]+1;
    for(int i=n-1;i>=0;i--){
        if(x[i]>m){
            if (i==0){
                cout<<"error_mono_next"<<endl;
            }else{
                x[i]=0;
                x[i-1]=x[i-1]+1;
            }
        }
    }
    return;
}

double CalculateMean(vector<double> & value)
{
    double sum = 0;
    int n=value.size();
    for(int i = 0; i < n; i++)
        sum += value[i];
    return (sum / n);
}

void getsubcell (const int m,const int k_order,vector<vector<int> > & v){
    if((m-1)>0){
        
        for (int i=0; i<=k_order;i++){
            
            vector<vector<int> > temp_v;
            getsubcell(m-1,k_order,temp_v);
            
            for(int j=0;j<temp_v.size();j++){
                temp_v[j].push_back(i);
            }
            for (int j=0;j<temp_v.size();j++){
                v.push_back(temp_v[j]);
            }
            
        }
    }
    else{
        for (int i=0; i<=k_order;i++){
            vector <int> v1;
            
            v1.push_back(i);
            v.push_back(v1);
        }
    }
    
}

void cell_generator (const int n,const int kmax, vector<vector<int> > & table ){
    //n: the number of variable
    //kmax: the maximum of interval's index(0,1,2,...,kmax)
    if(table.size() > 0){
        return;
    }
    
    vector<vector<int> > temp_table;
    getsubcell(n,kmax,temp_table);
    
    for (int j=0;j<temp_table.size();j++){
        vector<int> row;
        for (int k=temp_table[0].size();k>0;k--){
            row.push_back(temp_table[j][k-1]);
        }
        table.push_back(row);
    }
    temp_table.clear();
    
    
    
}

double log_mvnpdf(const vector<double> & mu,const  vector<double> & x,const vector <vector<double>>& inv_cov){
    //log Likelihood function : unnormalized value
    int n=mu.size();
    vector<double> dx(n);
    for(int i=0;i<n;i++){
        dx[i]=x[i]-mu[i];
    }
    //cout<<dx[0]<<" " <<dx[1]<<" " <<dx[2]<<" " <<dx[3]<<" " <<endl;
    vector<double> v(n,0.0);
    for(int i=0;i<inv_cov.size();i++){
        for(int j=0;j<n;j++){
            v[i]=v[i]+inv_cov[i][j]*dx[j];
        }
    }
    //cout<<inv_cov[0][0]<<" "<<inv_cov[0][1]<<" "<<inv_cov[0][2]<<" "<<inv_cov[0][3]<<endl;
//    cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<endl;
//    cout<<v[0]<<" "<<dx[1]<<" "<<dx[2]<<" "<<dx[3]<<endl;
    double res=0;
    for(int i=0;i<n;i++){
        res+=dx[i]*v[i];
    }
    

    res=-res/2.0;
    return res;
}

void computePosterior(const int numInt,const int sdim,const string model,const double dt,const int t,const vector<double> & cur_x,const vector<double> & next_x,const vector <vector<double>>& inv_cov,const vector <vector<double>> theta_list,map<string, FnPtr> & modelMap,vector<double> & logPrior,vector<double> & logPosterior){
    
    for(int i=0;i<theta_list.size();i++){
        
        vector<double> simNext_x = modelMap[model](cur_x,theta_list[i],dt);
        double logp=0;
//        for(int j=0;j<sdim;j++){
//            logp=logp-pow((next_x[j]-simNext_x[j]),2)/(2*cov[j][j]);
//        }


        //cout<<exp(log_mvnpdf(next_x,simNext_x,inv_cov))/sqrt(pow(2*pi,4)*0.000108736)<<endl;
        //if(t==3){cout<<exp(log_mvnpdf(next_x,simNext_x,inv_cov))<<endl;}
        logp=logp+log_mvnpdf(simNext_x,next_x,inv_cov);
        logPosterior[i]= logp+logPrior[i];

        
    }
    
    
    
    //normalization
    double normC=0.0,lognormC;
    for(int i=0;i<logPosterior.size();i++){
        if (isfinite(logPosterior[i])== false){ logPosterior[i]=-1e+12;}
        normC+=exp(logPosterior[i]);
    }
    lognormC=log(normC);
    if (isfinite(lognormC)== false){ lognormC=-1e+12;}
    for(int i=0;i<logPosterior.size();i++){
        logPosterior[i]-=lognormC;
    }

    
    
}
void bayesianPolyForm::UpAndLowGaussian_polyform(const int & j,const double & dy,const double &  maxdx,const double &  mindx,double & sum1_lo,double & sum1_up){
    //i:index of data time
    //j:index of observed states dim
    //double sum1=0.0;

    double temp_dy;


    //calculate the (upper/lower) bound of integration term

    if(dy > maxdx){
        temp_dy=dy-mindx;
        sum1_lo=exp(-pow(temp_dy,2)/(2.0 * cov[j]));
        temp_dy=dy-maxdx;
        sum1_up = exp(-pow(temp_dy,2)/(2.0 * cov[j]));
    }
    else if (dy < mindx){
        temp_dy=dy-maxdx;
        sum1_lo=exp(-pow(temp_dy,2)/(2.0 * cov[j]));
        temp_dy=dy-mindx;
        sum1_up = exp(-pow(temp_dy,2)/(2.0 * cov[j]));
    }
    else{
        sum1_up = 1.0;
        if(abs(dy-maxdx)>=abs(dy-mindx)){
            temp_dy=dy-maxdx;
        }
        else{
            temp_dy=dy-mindx;
        }
        sum1_lo=exp(-pow(temp_dy,2)/(2.0 * cov[j]));

    }

}

void bayesianPolyForm::CellTarget_pdf_polyform(const vector <double> & X, const vector <double> & delta_p,const double tol_uplo, vector <double> & logp)
{
    //tol_uplo: define the tolerance (p_up-p_lo/p_up) to refine the integration
    double p_lo=1.0;
    double p_up=1.0;
    double logp_lo,logp_up;

    int Np = X.size();
    //int d=(x0_initial.size()-Np);//dimension of state (not including parameters)
    //int indS=d+Np;
    int Ndeltak=9;//(Ndeltak+1)the initial number of point for integration to refine if (p_up-p_lo/p_up) > tol_uplo
    //int multiRef=2;//10//the multiply value for next refining iteration: Ndeltak=multiRef*previous(Ndeltak)
    //max number of iteration for integral
    //int maxIteration=2;//1000:No limit

    double cellMass=1.0;//the total weight of the cell: the total mass is equal for each cell. uniform prior:dp1*dp2..
    for (int i=0;i<delta_p.size();i++){
        cellMass = cellMass*(2.0*delta_p[i]);
    }


    //update the initial condition of ode (which considers parameters as states) with the current parameters
    int numSteps = tf/delta_t;
    //FitzhughNagumoPolyForm dpf(x0_initial[0], x0_initial[1], delta_t, numSteps, maxDegree);
    vector<MpfiWrapper> rangeVec;
//            MpfiWrapper range;
//            string varName1 = "alpha";
//            range.set(0.295, 0.305);
//            rangeVec.push_back(range);
//            string varName2 = "beta";
//            range.set(0.147, 0.153);
//            rangeVec.push_back(range);
//            vector<string> varNameVec;
//            varNameVec.push_back(varName1);
//            varNameVec.push_back(varName2);

    //vector<double> maxtheta_cur(Np), mintheta_cur(Np);
    for (int i=0; i<Np;i++){
        MpfiWrapper range;
        double mintheta_c= X[i]-delta_p[i];
        double maxtheta_c = X[i]+delta_p[i];
        range.set(mintheta_c, maxtheta_c);
        rangeVec.push_back(range);
    }
    computeNStepsAndSaveData(numSteps, maxDegree,varNameVec ,rangeVec);
    // get the flowpipes in x1 dimension
//    std::vector <double> result_upx1,result_lowx1;
//    dpf.getResultBound_x1(result_upx1,result_lowx1);

//    std::ofstream output("result_x.txt");
//    for (int i=0;i<globalSystem->result_upx1.size();i++)
//    {
//
//        output << globalSystem->result_upx1[i][0] << " "<<globalSystem->result_lowx1[i][0];
//        output << "\n";
//        //cout << globalSystem->result_upx1[i][0] << " "<<globalSystem->result_lowx1[i][0]<<endl;
//    }
//    output.close();
//    cout<<"size= "<<globalSystem->result_upx1.size()<<endl;

    int Nt = indSamp.size();//odesim1.data_y->size();
    for (int i=0;i<Nt;i++){
        //for (int i=1;i<2;i++){
        for (int j=0;j < obs_y[0].size();j++){
            int iS=indSamp[i];
            double dy = obs_y[i][j];
            //since result of flowpipe doesn't include state at t=0
            int iS_f=iS-1;
            double maxdx = globalSystem->result_upx1[iS_f][0];
            double mindx = globalSystem->result_lowx1[iS_f][0];
            double sum1_lo=0.0;
            double sum1_up=0.0;

            UpAndLowGaussian_polyform(j,dy, maxdx,mindx,sum1_lo,sum1_up);
//            if(1){
//                if(sum1_lo>sum1_up){
//                    cout<<"wrong  "<<maxdx<<" "<<mindx<<" "<<X[0]<<" "<<delta_p[0]<<X[1]<<" "<<delta_p[1]<<endl;
//                    cout<<"size= "<<result_upx1.size()<<endl;
//                    std::ofstream output("result_x.txt");
//                    for (int i1=0;i1<result_upx1.size();i1++)
//                    {
//
//                        output << result_upx1[i1] << " "<<result_lowx1[i1];
//                        output << "\n";
//                    }
//                    output.close();
//                }
//            }

            p_lo= p_lo * sum1_lo;
            p_up= p_up * sum1_up;
        }

    }
    rangeVec.clear();
//    result_upx1.clear();
//    result_lowx1.clear();


//    result.clear();
//    result.empty();
    p_lo= p_lo * cellMass;
    p_up=p_up * cellMass;
    cout<<p_lo<<" "<<p_up<<" " <<abs(p_up-p_lo)/abs(p_up)<<endl;
    //cout<<p_lo<<" "<<p_up<<" " <<p_up-p_lo<<endl;
    logp_up=log(p_up);
    logp_lo=log(p_lo);
    double dUpLo=abs(p_up-p_lo)/abs(p_up);
    double p_lo_old=p_lo;
    double p_up_old=p_up;
    double Ndeltak_temp=Ndeltak;

    int flag_p =0;
    if (isfinite(logp_up)== false){
        flag_p = 1;
        logp_lo=-1e+12;
        logp_up=-1e+12;
    }else{
        if (isfinite(logp_lo)== false){
            dUpLo=1.0;
        }
    }


    cout<<"finish"<<p_lo<<" "<<p_up<<endl;
    //    logp.push_back(log(p_lo));
    //    logp.push_back(log(p_up));
    logp.push_back(logp_lo);
    logp.push_back(logp_up);


}

void bayesianPolyForm::CellTarget_pdf_polyform_refine(const vector <double> & X, const vector <double> & delta_p,const double tol_uplo, vector <double> & logp){
    //tol_uplo: define the tolerance (p_up-p_lo/p_up) to refine the integration
    int Nt = indSamp.size();//odesim1.data_y->size();
    double logp_lo=logp[0],logp_up=logp[1];
    double p_lo=exp(logp_lo);
    double p_up=exp(logp_up);

    int Np = X.size();
    //int indS=d+Np;
    int Ndeltak=9;//(Ndeltak+1)the initial number of point for integration to refine if (p_up-p_lo/p_up) > tol_uplo
    int multiRef=2;//10//the multiply value for next refining iteration: Ndeltak=multiRef*previous(Ndeltak)
    //max number of iteration for integral
    int maxIteration=2;//1000:No limit

    double cellMass=1.0;//the total weight of the cell: the total mass is equal for each cell. uniform prior:dp1*dp2..
    for (int i=0;i<delta_p.size();i++){
        cellMass = cellMass*(2.0*delta_p[i]);
    }

    int numSteps = tf/delta_t;
    vector<MpfiWrapper> rangeVec;
//    vector<double> maxtheta_cur(Np), mintheta_cur(Np);
//    std::vector <double> result_upx1,result_lowx1;


    double dUpLo=abs(p_up-p_lo)/abs(p_up);
    double p_lo_old=p_lo;
    double p_up_old=p_up;
    double Ndeltak_temp=Ndeltak;

    int flag_p =0;
    if (isfinite(logp_up)!= false){
        if (isfinite(logp_lo)== false){
            dUpLo=1.0;
        }
    }

    int cout_I=1;
    if(!(dUpLo > tol_uplo && flag_p==0 && cout_I<maxIteration)){return;}
    while (dUpLo > tol_uplo && flag_p==0 && cout_I<maxIteration){
        cout_I=cout_I+1;
        //enumerate the number of monomials in n_var dimensions of degree up to max_order.
        int n_mono =  pow(Ndeltak_temp+1,Np);
        vector <int> f(Np);//initialize the order of monomial to be zeros
        p_lo=0.0;
        p_up=0.0;

        double cellMass_dm=cellMass/n_mono;
        vector <double> Celldelta_p, Cellmintheta;
        for (int i=0;i<delta_p.size();i++){
            double temp= delta_p[i]/double (Ndeltak_temp+1);//2.0*delta_p/(Ndeltak_temp+1)/2.0
            Celldelta_p.push_back(temp);
            Cellmintheta.push_back(X[i]-delta_p[i]);
        }

        for (int k=0;k<n_mono;k++){
            double p_lo_temp=1.0;
            double p_up_temp=1.0;

            //compute the flowpipe for each sub-cell in the current cell
            //update the initial condition of ode (which considers parameters as states) with the current parameters
            vector <double> temp_X(Np);
            //FitzhughNagumoPolyForm dpf1(x0_initial[0], x0_initial[1], temp_X[0], temp_X[1], Celldelta_p[0], Celldelta_p[1], delta_t, numSteps, maxDegree);
            for (int ik=0; ik<Np;ik++){
                MpfiWrapper range;
                double temp_X = Cellmintheta[ik]+Celldelta_p[ik]*(2*f[ik]+1);
                double mintheta_c = temp_X-Celldelta_p[ik];
                double maxtheta_c = temp_X+Celldelta_p[ik];
                range.set(mintheta_c, maxtheta_c);
                rangeVec.push_back(range);
            }
            // get the flowpipes in x1 dimension
            computeNStepsAndSaveData(numSteps, maxDegree,varNameVec ,rangeVec);


            //dpf.getResultBound_x1(result_upx1,result_lowx1);
//            cout<< "low "<<result_lowx1[49]<<endl;
//            cout<< "up "<<result_upx1[49]<<endl;

            for (int ii=0;ii<Nt;ii++){
                for (int j=0;j < obs_y[0].size();j++){
                    int iS=indSamp[ii];
                    double temp_dy = obs_y[ii][j];
                    //since result of flow* doesn't include state at t=0
                    int iS_f=iS-1;
                    double maxdx = globalSystem->result_upx1[iS_f][0];
                    double mindx = globalSystem->result_lowx1[iS_f][0];
                    double sum1_lo_temp=0.0;
                    double sum1_up_temp=0.0;
                    UpAndLowGaussian_polyform(j,temp_dy,maxdx,mindx,sum1_lo_temp,sum1_up_temp);
//                    if(sum1_lo_temp>sum1_up_temp){
//                        cout<<"wrong  "<<maxdx<<" "<<mindx<<" "<<endl;
//                        std::ofstream output("result_x.txt");
//                        for (int i=0;i<result_upx1.size();i++)
//                        {
//
//                            output << result_upx1[i] << " "<<result_lowx1[i];
//                            output << "\n";
//                        }
//                        output.close();
//                    }
//                    if(isnan(sum1_lo_temp)){
//                        cout<<"wrong  "<<maxdx<<" "<<mindx<<" "<<endl;
//                        cout<<"size_result  "<<result_upx1[49]<<endl;
//                    }
                    p_lo_temp= p_lo_temp * sum1_lo_temp;
                    p_up_temp= p_up_temp * sum1_up_temp;

                }

            }//Nt
            //cout<< "low_temp "<<p_lo_temp<<endl;
            p_lo_temp= p_lo_temp * cellMass_dm;
            p_up_temp= p_up_temp * cellMass_dm;
            p_lo=p_lo +p_lo_temp ;
            p_up=p_up +p_up_temp ;
//            cout<< "low "<<p_lo<<endl;
//            cout<< "up "<<p_up<<endl;
            if (k<n_mono-1){mono_next ( Ndeltak_temp, f );}
//            result_upx1.clear();
//            result_lowx1.clear();
            rangeVec.clear();
            //dpf1.~FitzhughNagumoPolyForm();
        }
        logp_up=log(p_up);
        logp_lo=log(p_lo);
        dUpLo=abs(p_up-p_lo)/abs(p_up);
        if (isfinite(logp_up)== false){
            logp_lo=-1e+12;
            logp_up=-1e+12;
            flag_p=1;
            //cout<<'flag'<<logp_up<<endl;
        }else{
            if (isfinite(logp_lo)== false){
                dUpLo=1.0;
            }
        }

        cout<<p_lo<<" "<<p_up<<" " <<abs(p_up-p_lo)/abs(p_up)<<endl;
        Ndeltak_temp=Ndeltak_temp*multiRef;

        //if(p_lo_old>p_lo||p_up_old<p_up){cout<<"wrong"<<endl;}

        p_lo_old=p_lo;
        p_up_old=p_up;


    }
    cout<<"finish"<<p_lo<<" "<<p_up<<endl;
    logp[0]=logp_lo;
    logp[1]=logp_up;

}


void bayesianPolyForm::allCell_polyform( const vector <double> & maxtheta,const vector <double> & mintheta, const double tol_uplo, vector <vector <double>> & theta_list, vector <double> & delta_p,vector <vector <double>> & thetaProb_list, vector <vector <double>> & res_deltap_list,vector <vector <double>> & res_theta_list,vector <vector <double>> & res_thetaProb_list){
    double threshold_integral = 3.0;
    vector <vector <int>> cell_table;
    for (int i=0;i<maxtheta.size();i++){
        double temp=(maxtheta[i]-mintheta[i])/numInt/2.0;
        delta_p.push_back(temp);
    }

    cell_generator (delta_p.size(),numInt-1,cell_table);



    vector <vector <double>> table_prob;
    vector <double> argMaxtheta;

    //double maximumProb = 0.0;
    vector <double> maximumProb(2,-100000000000.0);
    vector <double> normalcost(2, 0.0);


    //for (int i=0;i<5;i++){
    for (int i=0;i<cell_table.size();i++){
        //for (int i=0;i<1;i++){
        vector <double> theta_centra;
        for (int j=0; j<cell_table[0].size(); j++){
            double temp=mintheta[j]+delta_p[j]*(2*cell_table[i][j]+1);
            theta_centra.push_back(temp);
        }
//        cout<<"  "<<endl;
        vector <double> result_p;
        CellTarget_pdf_polyform(theta_centra,delta_p,tol_uplo, result_p);
        table_prob.push_back(result_p);



        if (result_p[0] > maximumProb[0]){
            argMaxtheta = theta_centra;
            maximumProb[0] = result_p[0];

        }
        if (result_p[1] > maximumProb[1]){maximumProb[1] = result_p[1];}
        //cout<<"p= "<<exp(table_prob[i][0])<<endl;

    }




    //find a set of maximum cells (consider the upper bound of prob)
    for (int i=0;i<cell_table.size();i++){
        vector <double> theta_centra;
        for (int j=0; j<cell_table[0].size(); j++){
            double temp=mintheta[j]+delta_p[j]*(2*cell_table[i][j]+1);
            theta_centra.push_back(temp);
        }

        //refine the integration by partitioning the cell with high prob into a set of sub-cells for computing an integral
        if (table_prob[i][1]+threshold_integral > maximumProb[1]) {
            cout<<"refine the integration"<<endl;
            CellTarget_pdf_polyform_refine(theta_centra, delta_p, tol_uplo, table_prob[i]);
        }


//        if (table_prob[i][1]+tor > maximumProb[1]){
//            //if (table_prob[i][0]+tor > 0.0){
//            theta_list.push_back(theta_centra);
//
//
//            thetaProb_list.push_back(table_prob[i]);
//
//        }
//        else{
//            res_theta_list.push_back(theta_centra);
//            res_deltap_list.push_back(delta_p);
//
//
//            res_thetaProb_list.push_back(table_prob[i]);
//
//        }
        res_theta_list.push_back(theta_centra);
        res_deltap_list.push_back(delta_p);
        res_thetaProb_list.push_back(table_prob[i]);

    }


}

void bayesianPolyForm::computeBayesian(const vector <double> & maxtheta,const vector <double> & mintheta){
    time_t tstart, tend;
    tstart = time(0);

    vector <double> delta_p;
    vector <vector <double>> theta_list;//[the n-th cell [p1 p2 ...]]
    vector <vector <double>> thetaProb_list;//[the n-th cell [lower, upper]]
    vector <vector <double>> res_deltap_list,res_theta_list,res_thetaProb_list;

    cout<<"iteration="<<0<<endl;
    double tol_uplo_temp=tol_uplo_lay[0];

    allCell_polyform(maxtheta,mintheta,tol_uplo_temp,theta_list, delta_p,thetaProb_list,res_deltap_list,res_theta_list,res_thetaProb_list);

//    for (int i_the=0;i_the<theta_list.size();i_the++){
//        cout<<theta_list[i_the][0]<<" "<<theta_list[i_the][1]<<" "<<"likelihood="<<thetaProb_list[i_the][0]<<", "<< thetaProb_list[i_the][1]<<endl;
//    }

    double maximumProb_res=-1e+12;
    for(int i=0;i<res_thetaProb_list.size();i++){
        if(res_thetaProb_list[i][0]>maximumProb_res){
            maximumProb_res=res_thetaProb_list[i][0];
        }
    }

    for(int i=0;i<res_thetaProb_list.size();i++){

        //prob multipy by a constant to avoid -inf computation
        res_thetaProb_list[i][0] = res_thetaProb_list[i][0] - maximumProb_res;
        res_thetaProb_list[i][1] = res_thetaProb_list[i][1] - maximumProb_res;
        //cout<<res_thetaProb_list[i][0]<< " "<<res_thetaProb_list[i][1]<<endl;
    }


    //calculate normalization constant
    double sumlo=0.0;
    double sumup=0.0;
    for(int i=0;i<res_thetaProb_list.size();i++){
        double temp_lo=exp(res_thetaProb_list[i][0]);
        double temp_up=exp(res_thetaProb_list[i][1]);
        if(isfinite(temp_up)==false){temp_up=0.0;}
        if(isfinite(temp_lo)==false){temp_lo=0.0;}
        sumlo=sumlo+temp_lo;
        sumup=sumup+temp_up;
    }

    //LP
    vector <vector <double>> res_thetaProb_listNorm=res_thetaProb_list;
    for(int i=0;i<res_thetaProb_list.size();i++){
        double w_up=exp(res_thetaProb_list[i][1]);
        double w_lo=exp(res_thetaProb_list[i][0]);
        if(isfinite(w_up)==false){w_up=0.0;}
        if(isfinite(w_lo)==false){w_lo=0.0;}

        double log_normWlo=log(w_lo/(w_lo-w_up+sumup));
        double log_normWup=log(w_up/(w_up-w_lo+sumlo));

        if (isfinite(log_normWlo)== false){log_normWlo=-1e+12;}
        if (isfinite(log_normWup)== false){log_normWup=-1e+12;}
        res_thetaProb_listNorm[i][0]=log_normWlo;
        res_thetaProb_listNorm[i][1]=log_normWup;
    }
    double sum1=0.0;
    double sum2=0.0;
    for (int i=0;i<res_thetaProb_listNorm.size();i++){
        sum1=sum1+exp(res_thetaProb_listNorm[i][0]);
        sum2=sum2+exp(res_thetaProb_listNorm[i][1]);
    }
    cout<<sum1<<" " <<sum2<<endl;



    tend = time(0);
    cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    //output result
    std::ofstream output_deltap("./outputs/deltap list.txt");
    for (int i=0;i<res_deltap_list.size();i++)
    {
        for (int j=0;j<res_deltap_list[0].size();j++)
        {
            output_deltap << res_deltap_list[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output_deltap << "\n";
    }
    output_deltap.close();

    std::ofstream output_theta("./outputs/centra list.txt");
    for (int i=0;i<res_theta_list.size();i++)
    {
        for (int j=0;j<res_theta_list[0].size();j++)
        {
            output_theta << res_theta_list[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output_theta << "\n";
    }
    output_theta.close();



    std::ofstream output_prob2("./outputs/prob list2.txt");
    for (int i=0;i<res_thetaProb_listNorm.size();i++)
    {
        //cout<<"likelihood="<<res_thetaProb_list[i][0]<<", "<< res_thetaProb_list[i][1]<<endl;
        for (int j=0;j<res_thetaProb_listNorm[0].size();j++)
        {
            output_prob2 << res_thetaProb_listNorm[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output_prob2 << "\n";
    }
    output_prob2.close();


}

void bayesianPolyForm::computeLikelihoodPolyFormBound(int maxDegree,const vector <double> & maxtheta,const vector <double> & mintheta) {
    MultivariatePoly logp(MpfiWrapper(0.0));
    StateAbstractionPtr st = fullSpace_pdf_polyform( logp);

    //computePolyIntegralAndSaveResults( maxtheta, mintheta, logp,st);
    computePolyIntegralAndSaveResults_Taylor( maxtheta, mintheta, logp,st);


}

void bayesianPolyForm::computePolyIntegralAndSaveResults(const vector <double> & maxtheta,const vector <double> & mintheta, const MultivariatePoly & logp, StateAbstractionPtr & st) {
    vector <double> delta_p;
    vector <vector <double>> thetaProb_list;//[the n-th cell [lower, upper]]
    vector <vector <double>> res_deltap_list,res_theta_list,res_thetaProb_list;
    //logp.prettyPrint(std::cout, std::map<int, string>());
    vector <vector <int>> cell_table;
    for (int i=0;i<maxtheta.size();i++){
        double temp=(maxtheta[i]-mintheta[i])/numInt/2.0;
        delta_p.push_back(temp);
    }

    double cellMass=1.0;//the total weight of the cell: the total mass is equal for each cell. uniform prior:dp1*dp2..
    for (int i=0;i<delta_p.size();i++){
        cellMass = cellMass*(2.0*delta_p[i]);
    }

    cell_generator (delta_p.size(),numInt-1,cell_table);

    vector <vector <double>> table_prob;
//    for (int i=0;i<5;i++){
    for (int i=0;i<cell_table.size();i++){
        //cout<<"iteration= "<<i<<endl;
        vector <double> theta_centra;
        for (int j=0; j<cell_table[0].size(); j++){
            double temp=mintheta[j]+delta_p[j]*(2*cell_table[i][j]+1);
            theta_centra.push_back(temp);
        }
//        cout<<"  "<<endl;
        vector <double> result_p;
        std::map<int, DistributionInfoPtr> Distrib_temp = globalSystem->getInitialMap();
        std::map<string, int> varIDs = globalSystem->getvarIDs();
        for(int j=0;j<theta_centra.size();j++){
            int varID = varIDs[varNameVec[j]];
            MpfiWrapper range(theta_centra[j],theta_centra[j]);
            /* Center the distribution? */
            MpfiWrapper offset = Distrib_temp[varID]->getOffset();
            range = range - offset;
            Distrib_temp[varID]->setRange(range);
        }
        std::map<int, MpfiWrapper> env_temp = ConvertToRangeMapForNoiseSymbols(Distrib_temp);
        MpfiWrapper logp_Mpfi= logp.evaluate(env_temp);
        result_p.push_back(logp_Mpfi.lower());//(logp_Mpfi.lower()+log(cellMass));
        result_p.push_back(logp_Mpfi.upper());//(logp_Mpfi.upper()+log(cellMass));
        table_prob.push_back(result_p);

        res_theta_list.push_back(theta_centra);
        res_deltap_list.push_back(delta_p);
    }
    res_thetaProb_list = table_prob;

    //find the maximum value for re-scaling
    double maximumProb_res=-1e+12;
    for(int i=0;i<res_thetaProb_list.size();i++){
        if(res_thetaProb_list[i][0]>maximumProb_res){
            maximumProb_res=res_thetaProb_list[i][0];
        }
    }

    for(int i=0;i<res_thetaProb_list.size();i++){

        //prob is multiplied by a constant to avoid -inf computation
        res_thetaProb_list[i][0] = res_thetaProb_list[i][0] - maximumProb_res;
        res_thetaProb_list[i][1] = res_thetaProb_list[i][1] - maximumProb_res;
        //cout<<res_thetaProb_list[i][0]<< " "<<res_thetaProb_list[i][1]<<endl;
    }

    vector <vector <double>> res_thetaProb_listNorm;
    normalizationOfInterval_discrete( res_thetaProb_list,  res_thetaProb_listNorm);

    //output result
    std::ofstream output_deltap("./outputs/deltap list.txt");
    for (int i=0;i<res_deltap_list.size();i++)
    {
        for (int j=0;j<res_deltap_list[0].size();j++)
        {
            output_deltap << res_deltap_list[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output_deltap << "\n";
    }
    output_deltap.close();

    std::ofstream output_theta("./outputs/centra list.txt");
    for (int i=0;i<res_theta_list.size();i++)
    {
        for (int j=0;j<res_theta_list[0].size();j++)
        {
            output_theta << res_theta_list[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output_theta << "\n";
    }
    output_theta.close();



    std::ofstream output_prob2("./outputs/prob list2.txt");
    for (int i=0;i<res_thetaProb_listNorm.size();i++)
    {
        //cout<<"likelihood="<<res_thetaProb_list[i][0]<<", "<< res_thetaProb_list[i][1]<<endl;
        for (int j=0;j<res_thetaProb_listNorm[0].size();j++)
        {
            output_prob2 << res_thetaProb_listNorm[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output_prob2 << "\n";
    }
    output_prob2.close();

}

void bayesianPolyForm::computePolyIntegralAndSaveResults_Taylor(const vector <double> & maxtheta,const vector <double> & mintheta, const MultivariatePoly & logp, StateAbstractionPtr & st) {
    //use Taylor series to approximate the polynomial exponential exp(p(x))
    logp.prettyPrint(std::cout, std::map<int, string>());
    MpfiWrapper logp_I = logp.getConstIntvl();//separate the interval from polynomial
    MultivariatePoly logp_poly = logp;
    logp_poly.setTerm( PowerProduct(),MpfiWrapper(0.0));//the polynomial without interval I
    //logp_poly.prettyPrint(std::cout, std::map<int, string>());
    MultivariatePoly prob = logp_poly.expPoly(st->getRangeMapForNoiseSymbols());
    prob = prob.truncate(maxDegree, st->getRangeMapForNoiseSymbols());
    MpfiWrapper exp_pI = exp(logp_I);
//    cout<<"I= "<<logp_I<<endl;
//    cout<<"expI= "<<exp_pI<<endl;
    //cout<<"prob = "<<endl;
    //prob.prettyPrint(std::cout, std::map<int, string>());
    prob.scaleAssign(exp_pI);
    //cout<<"after prob = "<<endl;
    //prob.prettyPrint(std::cout, std::map<int, string>());


//    std::map<int, DistributionInfoPtr> Distrib_temp = globalSystem->getInitialMap();
//    std::map<int, MpfiWrapper> env_temp = ConvertToRangeMapForNoiseSymbols(Distrib_temp);
//    MultivariatePoly prob = logp.expPoly(env_temp);
    //prob.prettyPrint(std::cout, std::map<int, string>());
    vector <double> delta_p;
    vector <vector <double>> thetaProb_list;//[the n-th cell [lower, upper]]
    vector <vector <double>> res_deltap_list,res_theta_list,res_thetaProb_list;

    vector <vector <int>> cell_table;
    for (int i=0;i<maxtheta.size();i++){
        double temp=(maxtheta[i]-mintheta[i])/numInt/2.0;
        delta_p.push_back(temp);
    }

    double cellMass=1.0;//the total weight of the cell: the total mass is equal for each cell. uniform prior:dp1*dp2..
    for (int i=0;i<delta_p.size();i++){
        cellMass = cellMass*(2.0*delta_p[i]);
    }

    cell_generator (delta_p.size(),numInt-1,cell_table);

    vector <vector <double>> table_prob;
//    for (int i=0;i<5;i++){
    for (int i=0;i<cell_table.size();i++){
        //cout<<"iteration= "<<i<<endl;
        vector <double> theta_centra;
        for (int j=0; j<cell_table[0].size(); j++){
            double temp=mintheta[j]+delta_p[j]*(2*cell_table[i][j]+1);
            theta_centra.push_back(temp);
        }
//        cout<<"  "<<endl;
        vector <double> result_p;
        std::map<int, DistributionInfoPtr> Distrib_temp = st->getNoiseSymbolInfoMap();//globalSystem->getInitialMap();
        std::map<string, int> varIDs = globalSystem->getvarIDs();
        for(int j=0;j<theta_centra.size();j++){
            int varID = varIDs[varNameVec[j]];
            //MpfiWrapper range(theta_centra[j]-delta_p[j],theta_centra[j]+delta_p[j]);
            MpfiWrapper range(theta_centra[j],theta_centra[j]);
            /* Center the distribution? */
            MpfiWrapper offset = Distrib_temp[varID]->getOffset();
            range = range - offset;
            Distrib_temp[varID]->setRange(range);

            //move to center zero
//            MpfiWrapper offset_tmp = offset + median(range);
//            Distrib_temp[varID]->setOffset(offset_tmp);
//            range = range - median(range);
//            Distrib_temp[varID]->setRange(range);

            //cout<< j<<" range = " <<range<<endl;
        }
        std::map<int, MpfiWrapper> env_temp = ConvertToRangeMapForNoiseSymbols(Distrib_temp);

//        MpfiWrapper logp_I = logp.getConstIntvl();//separate the interval from polynomial
//        MultivariatePoly logp_poly = logp;
//        logp_poly.setTerm( PowerProduct(),MpfiWrapper(0.0));//the polynomial without interval I
//        //logp_poly.prettyPrint(std::cout, std::map<int, string>());
//        MultivariatePoly prob = logp_poly.expPoly(env_temp);

        MpfiWrapper prob_Mpfi= prob.evaluate(env_temp);//logp.evaluate(env_temp);

        result_p.push_back(log(prob_Mpfi.lower()));//(logp_Mpfi.lower()+log(cellMass));
        result_p.push_back(log(prob_Mpfi.upper()));//(logp_Mpfi.upper()+log(cellMass));
        table_prob.push_back(result_p);
        cout<<prob_Mpfi.lower()<<", "<<prob_Mpfi.upper()<<endl;
        res_theta_list.push_back(theta_centra);
        res_deltap_list.push_back(delta_p);
    }
    res_thetaProb_list = table_prob;

    //find the maximum value for re-scaling
    double maximumProb_res=-1e+12;
    for(int i=0;i<res_thetaProb_list.size();i++){
        if(res_thetaProb_list[i][0]>maximumProb_res){
            maximumProb_res=res_thetaProb_list[i][0];
        }
    }

    for(int i=0;i<res_thetaProb_list.size();i++){

        //prob is multiplied by a constant to avoid -inf computation
        res_thetaProb_list[i][0] = res_thetaProb_list[i][0] - maximumProb_res;
        res_thetaProb_list[i][1] = res_thetaProb_list[i][1] - maximumProb_res;
        //cout<<res_thetaProb_list[i][0]<< " "<<res_thetaProb_list[i][1]<<endl;
    }

    vector <vector <double>> res_thetaProb_listNorm;
    normalizationOfInterval_discrete( res_thetaProb_list,  res_thetaProb_listNorm);

    //output result
    std::ofstream output_deltap("./outputs/deltap list.txt");
    for (int i=0;i<res_deltap_list.size();i++)
    {
        for (int j=0;j<res_deltap_list[0].size();j++)
        {
            output_deltap << res_deltap_list[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output_deltap << "\n";
    }
    output_deltap.close();

    std::ofstream output_theta("./outputs/centra list.txt");
    for (int i=0;i<res_theta_list.size();i++)
    {
        for (int j=0;j<res_theta_list[0].size();j++)
        {
            output_theta << res_theta_list[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output_theta << "\n";
    }
    output_theta.close();



    std::ofstream output_prob2("./outputs/prob list2.txt");
    for (int i=0;i<res_thetaProb_listNorm.size();i++)
    {
        //cout<<"likelihood="<<res_thetaProb_list[i][0]<<", "<< res_thetaProb_list[i][1]<<endl;
        for (int j=0;j<res_thetaProb_listNorm[0].size();j++)
        {
            output_prob2 << res_thetaProb_listNorm[i][j] << " "; // behaves like cout - cout is also a stream
        }
        output_prob2 << "\n";
    }
    output_prob2.close();

}


void bayesianPolyForm::test(int n){
//    for (int i = 0 ; i < 2; ++i){
//        std::cout << "After " << i << " Steps" << std::endl;
//        globalSystem -> prettyPrintStateAbstraction(std::cout, st);
//        computeNSteps_inference(n, maxDegree, st);
//        //evaluateLogLikelihood_AtCur( st,maxDegree,0.2,);
//    }
    MultivariatePoly logp(MpfiWrapper(0.0));
    StateAbstractionPtr st = fullSpace_pdf_polyform( logp);
    //std::map<int, DistributionInfoPtr> Distrib_temp = globalSystem->getInitialMap();
    //std::map<int, MpfiWrapper> env_temp = ConvertToRangeMapForNoiseSymbols(Distrib_temp);
    //MultivariatePoly result = logp.expPoly(env_temp);
    logp.prettyPrint(std::cout, std::map<int, string>());
    //result.prettyPrint(std::cout, std::map<int, string>());

}
std::map<int, MpfiWrapper> bayesianPolyForm::ConvertToRangeMapForNoiseSymbols(std::map<int, DistributionInfoPtr> & Distrib){
    std::map<int, MpfiWrapper> retMap;
    for (auto p: Distrib){
        retMap.insert(make_pair(p.first, p.second->getRange()));

    }
    return retMap;
}

StateAbstractionPtr bayesianPolyForm::fullSpace_pdf_polyform( MultivariatePoly & logp)
{
//    logp.setConst(0.0);

//    int Np = X.size();
    //int d=(x0_initial.size()-Np);//dimension of state (not including parameters)
    //int indS=d+Np;


    //update the initial condition of ode (which considers parameters as states) with the current parameters
    int numSteps = Tdel/delta_t;//tf/delta_t;

    StateAbstractionPtr st = initialize_globalSystem( maxDegree);
    int Nt = indSamp.size();//odesim1.data_y->size();
    for (int i=0;i<Nt;i++){
        //for (int i=1;i<2;i++){
        for (int j=0;j < obs_y[0].size();j++){
            //int iS=indSamp[i];
            double dy = obs_y[i][j];
            //since result of flowpipe doesn't include state at t=0
            //int iS_f=iS-1;
            computeNSteps_inference(numSteps, maxDegree, st);
            evaluateLogLikelihood_AtCur( st,maxDegree,  dy, logp);
        }

    }

    //globalSystem -> prettyPrintStateAbstraction(std::cout, st);
    //logp.prettyPrint(cout, std::map<int, string> ());
    return st;

}

void bayesianPolyForm::evaluateLogLikelihood_AtCur( StateAbstractionPtr st,int maxDegree, double dy, MultivariatePoly & logp) {
    //std::map<int, string> name_env;
    std::vector<Query> queries;
    queries = globalSystem -> getQueries();
    for (auto q: queries){
        ExprPtr e = q.getExpr();
        MultivariatePoly p = e -> evaluate(st);
        p.scaleAndAddAssign(-1.0 , MultivariatePoly(MpfiWrapper(dy)));
        MultivariatePoly p2=p.squarePoly();
        p2.scaleAssign(MpfiWrapper(1.0/(2.0 * cov[0])));
        //if(maxDegree > 0 ) {
        MultivariatePoly p2_trunc = p2.truncate(maxDegree, st->getRangeMapForNoiseSymbols());
        p2_trunc.centerAssign(st->getRangeMapForNoiseSymbols());
           // p2_trunc.prettyPrint(what, name_env);
        //}
        //logp.scaleAndAddAssign(-1.0, p2);
        logp.scaleAndAddAssign(-1.0, p2_trunc);

    }


}

void bayesianPolyForm::normalizationOfInterval_discrete( vector <vector <double>> const & res_thetaProb_list, vector <vector <double>> & res_thetaProb_listNorm){
    //calculate normalization constant
    double sumlo=0.0;
    double sumup=0.0;
    for(int i=0;i<res_thetaProb_list.size();i++){
        double temp_lo=exp(res_thetaProb_list[i][0]);
        double temp_up=exp(res_thetaProb_list[i][1]);
        if(isfinite(temp_up)==false){temp_up=0.0;}
        if(isfinite(temp_lo)==false){temp_lo=0.0;}
        sumlo=sumlo+temp_lo;
        sumup=sumup+temp_up;
    }

    //LP
    res_thetaProb_listNorm=res_thetaProb_list;
    for(int i=0;i<res_thetaProb_list.size();i++){
        double w_up=exp(res_thetaProb_list[i][1]);
        double w_lo=exp(res_thetaProb_list[i][0]);
        if(isfinite(w_up)==false){w_up=0.0;}
        if(isfinite(w_lo)==false){w_lo=0.0;}

        double log_normWlo=log(w_lo/(w_lo-w_up+sumup));
        double log_normWup=log(w_up/(w_up-w_lo+sumlo));

        if (isfinite(log_normWlo)== false){log_normWlo=-1e+12;}
        if (isfinite(log_normWup)== false){log_normWup=-1e+12;}
        res_thetaProb_listNorm[i][0]=log_normWlo;
        res_thetaProb_listNorm[i][1]=log_normWup;
    }
    double sum1=0.0;
    double sum2=0.0;
    for (int i=0;i<res_thetaProb_listNorm.size();i++){
        sum1=sum1+exp(res_thetaProb_listNorm[i][0]);
        sum2=sum2+exp(res_thetaProb_listNorm[i][1]);
    }
    cout<<sum1<<" " <<sum2<<endl;
}