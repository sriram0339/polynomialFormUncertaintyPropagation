//
// Created by Sriram Sankaranarayanan on 2/14/20.
//
#include "ModelParser.hh"

StochasticSystem * globalSystem;
extern "C" FILE *yyin;
extern int yyparse();

void parserMain(const char * fileName){
    globalSystem = new StochasticSystem();
    FILE * f = fopen(fileName, "r");
    if (!f) {
        std::cerr << "FATAL: failed to open file -- " << fileName << std::endl;
        exit(1);
    }

    yyin = f;
    yyparse();
    fclose(f);
    return;
}



