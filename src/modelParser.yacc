%{
	#include "ModelParser.hh"
	#include "modelParser.yy.hh"
	extern int yyerror(const char *);
	extern int yyerror(std::string);
	extern int yylex();
	extern int yyparse();
	bool err;

	int lineNum = 1;
	void parseError(const char *str, int lnum)
    	{
    		std::cerr << "Error @line " << lineNum << ":" << std::string(str) << std::endl;
    		exit(1);
    	}

    using namespace PolynomialForms;

%}

%union
{
    string * identifier;
    double dblVal;
    int intval;
    VariableDistributionInfo * dPtr;
    Expr * expr;
}

%token <dblVal> NUM
%token <identifier> IDENT
%token ASSIGN
%token INIT
%token BELONGSTO
%token MSIN
%token MCOS
%token MTAN
%token MSQRT
%token MEXP
%token VAR
%token DISTURB
%token END
%token COLON
%token VISUALIZE
%token QUERY
%token UNIFORM
%token EXPECT
%token TRUNCGAUSSIAN
%token GEQ
%token LEQ

%type<dPtr> distributionSpec
%type<expr> expr
%type<intval> identifier

%right '^'
%left '+' '-'
%left '*' '/'
%nonassoc uminus
%left ASSIGN


%start model

%%

model: varDeclarations modelEquations queries END;

varDeclarations: singleVarDecl
| varDeclarations singleVarDecl;



singleVarDecl: VAR IDENT  INIT distributionSpec
{

    DistributionInfoPtr dPtr = std::shared_ptr<VariableDistributionInfo>($4);
    globalSystem -> addVar(*$2, dPtr);
    delete($2);
};

distributionSpec: UNIFORM '(' NUM ',' NUM ')' {
    $$ = new UniformDistributionInfo(MpfiWrapper($3, $5));
}
| TRUNCGAUSSIAN '(' NUM ',' NUM ',' NUM ',' NUM ')' {
    $$ = new TruncNormalDistributionInfo(MpfiWrapper($7, $9), $3, $5);
}
;


modelEquations: modelEquations equation
| equation;

equation: identifier ASSIGN expr {
    globalSystem -> addUpdate($1, ExprPtr($3));
};

expr: expr '+' expr {
    $$ = new Plus( ExprPtr($1), 1.0, ExprPtr($3), 1.0);
}

| expr '-' expr {
    $$ = new Plus( ExprPtr($1), 1.0, ExprPtr($3), -1.0);
}

| expr '/' expr {
    $$ = new Div(ExprPtr($1), ExprPtr($3));
}

| expr '*' expr {
    $$ = new Star(ExprPtr($1), ExprPtr($3));
}

| expr '^' NUM {
    $$ = new Pow(ExprPtr($1), (int) $3);
}

| MSIN '(' expr ')' {
    $$ = new Trig(SIN_TYPE, ExprPtr($3));
}

| MCOS '(' expr ')' {
    $$ = new Trig(COS_TYPE, ExprPtr($3));
}


| '(' expr ')' { $$ = $2; }
| '-' expr %prec uminus {
    Plus * p  = new Plus();
    p -> addFactor(ExprPtr($2), -1.0);
    $$ = p;
}

| NUM {
    $$ = new Const($1);

 }
| identifier {
    $$ = new Var($1);
 }
| distributionSpec {
    $$ = new Distrib(DistributionInfoPtr($1));
};

identifier : IDENT {
    $$ = globalSystem -> getVarIDFromName(*$1);
    if ($$ < 0){
        std::cerr << "FATAL: could not find " << *$1 << " at line number: " << lineNum << std::endl;
        yyerror("Variable not found");
    }
    delete($1);
 };

queries: queries singleQuery
| singleQuery;

singleQuery: QUERY IDENT expr GEQ expr  {
    ExprPtr queryExpr = std::make_shared<Plus>(ExprPtr($3), 1.0, ExprPtr($5), -1.0);
    Query q = probabilityQuery(*$2, queryExpr);
    delete($2);
    globalSystem -> addQuery(q);
}
| QUERY IDENT expr LEQ expr {
    ExprPtr queryExpr = std::make_shared<Plus>(ExprPtr($3), -1.0, ExprPtr($5), 1.0);
    Query q = probabilityQuery(*$2, queryExpr);
    delete($2);
    globalSystem -> addQuery(q);
}
| QUERY IDENT EXPECT expr {
    ExprPtr queryExpr($4);
    Query q = expectationQuery(*$2, queryExpr);
    globalSystem -> addQuery(q);
    delete($2);
};


%%

extern int yyerror(char const * err){
    parseError(err, lineNum);
    exit(1);
    return -1;
}