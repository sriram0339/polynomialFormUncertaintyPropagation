//
// Created by Sriram Sankaranarayanan on 2/10/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_EXPRAST_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_EXPRAST_HH

#include <string>
#include <vector>
#include <map>
#include "MpfiWrapper.hh"
#include "DistributionInfo.hh"
#include "MultivariatePoly.hh"

namespace PolynomialForms {
    typedef enum {VAR_TYPE, CONST_TYPE, PLUS_TYPE, MULT_TYPE, SIN_TYPE, COS_TYPE, DISTRIB_TYPE, POW_TYPE} expr_type_t;

    class ExprVisitor;
    typedef std::shared_ptr<ExprVisitor> ExprVisitorPtr;
    
    class Expr {
    protected:
        expr_type_t eType;
    public:
        explicit Expr(expr_type_t ty_): eType(ty_) {};
        [[nodiscard]] expr_type_t getType() const { return eType; };
        virtual void visit(ExprVisitorPtr ev) const = 0;
        virtual MultivariatePoly evaluate(std::map<int, MultivariatePoly> const & env,
                                          std::map<int, DistributionInfoPtr>  & randomVarInstancesMap)  const = 0;

        virtual ~Expr() = default;
    };

    typedef shared_ptr<Expr> ExprPtr;

    class Const: public Expr {
    protected:
        MpfiWrapper c;
    public:
        Const(MpfiWrapper c_): Expr(CONST_TYPE), c(c_) {};
        virtual void visit(ExprVisitorPtr ev ) const;
        virtual MultivariatePoly evaluate(std::map<int, MultivariatePoly> const & env,
                                          std::map<int, DistributionInfoPtr>  & randomVarInstancesMap) const {
            return MultivariatePoly(c);
        }
    };

    class Var: public Expr {
    protected:
       int varID;
    public:
        Var(int id_): Expr(VAR_TYPE), varID(id_) {};
        virtual void visit(ExprVisitorPtr ev ) const;
        virtual MultivariatePoly evaluate(std::map<int, MultivariatePoly> const & env,
                                          std::map<int, DistributionInfoPtr>  & randomVarInstancesMap) const {
            auto it = env.find(varID);
            assert(it != env.end());
            return it -> second;
        }
    };

    class Plus: public Expr {
    protected:
        std::vector<ExprPtr> subExprs;
        std::vector<MpfiWrapper> scaleFactors;
    public:
        Plus(ExprPtr e1, MpfiWrapper c1, ExprPtr e2, MpfiWrapper c2):
        Expr(PLUS_TYPE),
        subExprs({e1, e2}),
        scaleFactors({c1, c2}){};

        Plus(): Expr(PLUS_TYPE){};

        void addFactor(ExprPtr e, MpfiWrapper scale = MpfiWrapper(1.0)) {
            subExprs.push_back(e);
            scaleFactors.push_back(scale);
        }

        virtual void visit(ExprVisitorPtr ev ) const ;
        virtual MultivariatePoly evaluate(std::map<int, MultivariatePoly> const & env,
                                          std::map<int, DistributionInfoPtr>  & randomVarInstancesMap) const{
          int i =0;
          MultivariatePoly retPoly;
          for (i = 0; i < subExprs.size(); ++i ){
              MultivariatePoly tmp = subExprs[i] -> evaluate(env, randomVarInstancesMap);
              retPoly.scaleAndAddAssign(scaleFactors[i], tmp);
          }
          return retPoly;
        }
    };

    class Star: public Expr {
    protected:
        std::vector<ExprPtr> subExprs;
    public:
        Star(ExprPtr e1, ExprPtr e2): Expr(MULT_TYPE), subExprs({e1, e2}) {};


        virtual void visit(ExprVisitorPtr ev ) const;
        virtual MultivariatePoly evaluate(std::map<int, MultivariatePoly> const & env,
                                          std::map<int, DistributionInfoPtr>  & randomVarInstancesMap) const{
            MultivariatePoly retPoly(1.0);
            int i = 0;
            for (i = 0; i < subExprs.size(); ++i ) {
                MultivariatePoly tmp = subExprs[i] -> evaluate(env, randomVarInstancesMap);
                retPoly = retPoly.multiply(tmp);
            }
            return retPoly;
        }
    };

    class Pow: public Expr {
    protected:
        ExprPtr subExpr;
        int pow;
    public:
        Pow(ExprPtr e, int j): Expr(POW_TYPE), subExpr(e), pow(j){
            assert(pow >= 1);
        };

        virtual void visit(ExprVisitorPtr ev ) const ;
        virtual MultivariatePoly evaluate(std::map<int, MultivariatePoly> const & env,
                                          std::map<int, DistributionInfoPtr>  & randomVarInstancesMap) const {
            MultivariatePoly tmp = subExpr -> evaluate(env, randomVarInstancesMap);
            return tmp.powPoly(pow);
        }
    };

    class Trig: public Expr {
    protected:
        ExprPtr subExpr;
    public:
        Trig(expr_type_t eType, ExprPtr e): Expr(eType), subExpr(e) {
            assert(eType == SIN_TYPE || eType == COS_TYPE);
        }
        virtual void visit(ExprVisitorPtr ev ) const ;
        virtual MultivariatePoly evaluate(std::map<int, MultivariatePoly> const & env,
                std::map<int, DistributionInfoPtr>  & randomVarInstancesMap) const {
            MultivariatePoly tmp = subExpr -> evaluate(env, randomVarInstancesMap);
            std::map<int, MpfiWrapper> rangeMap;
            for (auto p: randomVarInstancesMap){
                rangeMap.insert(std::make_pair(p.first, p.second -> getRange()));
            }
            // to do
            switch(eType){
                case SIN_TYPE:
                    return tmp.sine(rangeMap);
                    break;
                case COS_TYPE:
                    return tmp.cosine(rangeMap);
                    break;
                default:
                    assert(false);
            }
        }
    };

    class Distrib: public Expr {
    protected:
        DistributionInfoPtr  dPtr;
    public:
        Distrib(DistributionInfoPtr dPtr_): Expr(DISTRIB_TYPE), dPtr(dPtr_) {};
        virtual void visit(ExprVisitorPtr ev) const;
        virtual MultivariatePoly evaluate(std::map<int, MultivariatePoly> const & env,
                                          std::map<int, DistributionInfoPtr>  & randomVarInstancesMap) const {
            // Make a new instance
            int varId = randomVarInstancesMap.size();
            MultivariatePoly retPoly(1.0, varId);
            randomVarInstancesMap.insert(std::make_pair(varId, dPtr));
            return retPoly;
        }
    };


    class ExprVisitor{
    public:
        virtual void visitDistrib(DistributionInfoPtr dPtr){}
        virtual void visitTrig(expr_type_t eType, ExprPtr subExpr){}
        virtual void visitStar(std::vector<ExprPtr> const & subExprs){}
        virtual void visitPlus(std::vector<ExprPtr> const & subExprs, std::vector<MpfiWrapper> const & scaleFactors){}
        virtual void visitVar(int varID){}
        virtual void visitConst(MpfiWrapper const & c){}
        virtual void visitPow(ExprPtr subExpr, int pow){}
    };





};


#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_EXPRAST_HH
