//
// Created by Sriram Sankaranarayanan on 2/10/20.
//

#include "ExprAST.hh"

namespace PolynomialForms {
    void Const::visit(ExprVisitorPtr ev) const {
        ev->visitConst(c);
    }


    void Var::visit(ExprVisitorPtr ev) const {
        ev -> visitVar(varID);
    }


    void Plus::visit(ExprVisitorPtr ev) const {
        ev -> visitPlus(subExprs, scaleFactors);
    }


    void Star::visit(ExprVisitorPtr ev) const {
        ev -> visitStar(subExprs);
    }

    void Trig::visit(ExprVisitorPtr ev) const {
        ev -> visitTrig(eType, subExpr);
    }

    void Distrib::visit(ExprVisitorPtr ev) const {
        ev -> visitDistrib(dPtr);
    }

    void Pow::visit(ExprVisitorPtr ev) const {
        ev -> visitPow(subExpr, pow);
    }

    void Div::visit(ExprVisitorPtr ev) const {
        ev -> visitDiv(numer, denom);
    }
}


