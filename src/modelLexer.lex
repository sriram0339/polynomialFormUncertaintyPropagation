%{
/*---
 Lexcical Analysis
---*/
#include <cstdlib>
#include "ModelParser.hh"
#include "modelParser.yy.hh"
%}

delim [ \t\r]
line [\n]
whitespace {delim}+
digit [0-9]
letter [a-zA-Z_]
ident {letter}({letter}|{digit}|".")*
number ("-"?)(({digit}+)|({digit}*"."{digit}*)|({digit}+"e"(({digit}+)|("-"{digit}+)))|({digit}*"."{digit}*"e"(({digit}+)|("-"{digit}+)))|({digit}*"."{digit}*"e"("+"({digit}+)|("-"{digit}+))))
stringLiteral \"(\\.|[^"\\])*\"

%%
"+" {return '+';}
"*" {return '*';}
"-" {return '-';}
"," {return ',';}
";" {return ';';}
":" {return ':';}
"(" {return '(';}
")" {return ')';}
"{" {return '{';}
"}" {return '}';}
"[" {return '[';}
"]" {return ']';}
":=" {return ASSIGN;}
"^" {return '^';}
"/" {return '/';}
">=" {return GEQ; }
"<=" {return LEQ; }
"init" {return INIT;}
"end" {return END;}
"exp" {return MEXP;}
"sin" {return MSIN;}
"cos" {return MCOS;}
"sqrt" {return MSQRT;}
"var" {return VAR;}
"query" {return QUERY; }
"expectation" {return EXPECT; }
"uniform" {return UNIFORM; }
"truncGaussian" {return TRUNCGAUSSIAN; }
"observed" {return OBSERV; }


{number} { yylval.dblVal = atof( (char *)yytext ); return NUM; }

{ident}	{ yylval.identifier = new std::string(yytext); return IDENT; }

{stringLiteral} {
    char * tmpString = new char[strlen(yytext)+1];
    strcpy(tmpString, yytext+1);
    tmpString[strlen(yytext)-2] = '\0';
    std::cout << "String literal: " << tmpString << std::endl;
    yylval.identifier = new std::string(tmpString);
    delete[](tmpString);
    return IDENT;
   }

{whitespace}

{line} { lineNum++; }

"#" {	/* Comment line */
	int c;
	c = yyinput();
	while(c!='\n' && c!=EOF)
	{
		c = yyinput();
	}

	if(c == '\n')
	{
		++lineNum;
	}
}

.  { printf("Unknown Character in line %d : %s -- Ignored\n", lineNum, yytext); }

%%


int yywrap()
{
	return 1;
}

