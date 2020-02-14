//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_QUERIES_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_QUERIES_HH

#include <string>
#include "ExprAST.hh"

namespace PolynomialForms {

    class Query {
    protected:
        std::string id;
        ExprPtr queryExpr; // we seek its expectation
    public:
        Query(std::string qid, ExprPtr e_) : id(qid), queryExpr(e_) {}
    };

    class ProbabilityQuery: public Query {
    protected:
    public:
        ProbabilityQuery(std::string qid, ExprPtr e_): Query(qid, e_){};
    };

    class ExpectationQuery: public Query {
    protected:

    public:
        ExpectationQuery(std::string qid, ExprPtr e_): Query(qid, e_){};
    };

};
#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_QUERIES_HH
