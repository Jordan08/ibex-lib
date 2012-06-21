%{
//============================================================================
//                                  I B E X                                   
// File        : Yacc/Bison input for Ibex parser
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jun 12, 2012
// Last Update : Jun 12, 2012
//==========================

#include <math.h>
#include <vector>

#include "ibex_ParserExpr.h"
#include "ibex_ParserNumConstraint.h"
#include "ibex_ParserGenerator.h"
#include "ibex_Array.h"

using namespace std;

extern int ibexlex(void);
extern char* ibextext;
extern int ibex_lineno;

void ibexerror (const std::string& msg) {
  throw ibex::SyntaxErrorException(msg, ibextext, ibex_lineno);
}
   
namespace ibex { 
namespace parser {

P_Source source;
P_result result;

stack<Scope> scopes;

void begin() {
  ibex_lineno=-1;
  if (!setlocale(LC_NUMERIC, "C")) // to accept the dot (instead of the french coma) with numeric numbers
    ibexerror("platform does not support \"C\" locale");
  
  ibex_lineno=1;

  /* there may be some pending scopes (if the previous call to the parser failed).
   */
  while (!scopes.empty()) scopes.pop(); 
  
  scopes.push(Scope()); // a fresh new scope!
}

int _2int(const ExprNode& expr) {
	return ConstantEvaluator(scopes.top()).eval_integer(expr);
}

double _2dbl(const ExprNode& expr) {
	return ConstantEvaluator(scopes.top()).eval_double(expr);
}

void end() {
	Generator(source,result);
}

} // end namespace
} // end namespace


using namespace ibex;
using namespace parser;

%}	

%union{
  char*     str;
  int       itg;
  double    real;
  Interval* itv;
  vector<const char*>*                   strvec;

  ibex::Dim*                             dim;
  const ibex::parser::P_ExprSymbol*      sbl;

  vector<ibex::parser::P_NumConstraint*> _ctrblklst;
  ibex::parser::P_NumConstraint*         _ctrblk;
  const ibex::NumConstraint*             _ctr;

  vector<ibex::parser::P_ConstantExpr>* cstvec;
  ibex::parser::P_ConstantExpr*         cst;
  const ibex::ExprNode*                 _expr;
  vector<const ibex::ExprNode*>*        exprvec;

}

%token <str> TK_NEW_SYMBOL
%token <str> TK_TMP_SYMBOL
%token <str> TK_INPUT_SYMBOL
%token <str> TK_RETURN_SYMBOL
%token <str> TK_ITERATOR
%token <str> TK_CTR_LIST
%token <str> TK_CTR
%token <str> TK_FUNC

%token <str> TK_STRING
%token <itg> TK_INTEGER
%token <real> TK_FLOAT
%token <_bool> TK_BOOL
%token <str> TK_CONSTANT

%token TK_PARAM TK_CONST TK_VARS TK_FUNCTION 
%token TK_MIN TK_MAX TK_INF TK_MID TK_SUP TK_SIGN TK_ABS
%token TK_SQRT TK_EXPO  TK_LOG 
%token TK_COS  TK_SIN   TK_TAN  TK_ARCCOS  TK_ARCSIN  TK_ARCTAN TK_ARCTAN2
%token TK_COSH TK_SINH  TK_TANH TK_ARCCOSH TK_ARCSINH TK_ARCTANH
%token TK_LEQ  TK_GEQ   TK_EQU  TK_ASSIGN

%token TK_BEGIN TK_END TK_FOR TK_FROM TK_TO

%token TK_CTRS TK_CONSTRAINT TK_CONSTRAINT_LIST

%nonassoc TK_EQU

%nonassoc '<' TK_LEQ '>' TK_GEQ    /* on interdit donc a < b < c */
%left '+' '-' TK_UNION
%left '*' '/' TK_INTERSEC
%left '^'


/* --------------------------- Nonterminals types -------------------------- */

%type<sbl>        identifier
%type<dim>        dimension

/* functions */
%type<strvec>     fnc_inpt_list
%type<str>        fnc_input
%type<ass>        fnc_assign

/* symbols */
%type<sbl>        decl_sbl

/* constraints */
%type<_ctrblklst> decl_ctr_list
%type<_ctrblklst> ctr_blk_list
%type<_ctrblklst> ctr_blk_list_ // ctr_blk_list without ending semicolon
%type<_ctrblk>    ctr_blk
%type<_ctrblk>    ctr_loop
%type<_ctr>       ctr 
 
%type<itv>        interval

/* expressions */
%type<_expr>      expr
%type<exprvec>    expr_list

%%


program       :                                      { begin(); }
              decl_opt_cst 
              TK_VARS                   
	          decl_var_list ';'   
              decl_opt_par
              decl_fnc_list
	          decl_ctr_list                          { end(); }
              ;

/**********************************************************************************************************************/
/*                                                SYMBOLS                                                         */
/**********************************************************************************************************************/

decl_opt_cst  : 
              | TK_CONST decl_cst_list ';'
	          ;

decl_cst_list : decl_cst                                    
              | decl_cst_list ';' decl_cst
              ;

decl_cst      : TK_NEW_SYMBOL TK_EQU const_expr      { scopes.top().add_cst($1, *$3); free($1); delete $3; }
              | TK_NEW_SYMBOL TK_IN const_expr       { scopes.top().add_cst($1, *$3); free($1); delete $3; }
              ;

decl_opt_par  : 
              | TK_PARAM decl_par_list ';'
              ; 
 
decl_par_list : decl_par                                    
              | decl_par_list ';' decl_par
              | decl_par_list ',' decl_par
              ;

decl_par      : '!' decl_sbl                         { source.eprs.push_back($1); }
	          |     decl_sbl                         { source.sybs.push_back($1); }
              ;

decl_var_list : decl_var                             
              | decl_var_list ';' decl_sbl           
              | decl_var_list ',' decl_sbl           
              ;

decl_var      : decl_sbl                             { source.vars.push_back($1); }
              ;
              
decl_sbl      : TK_NEW_SYMBOL dimension              { $$ = new P_ExprSymbol($1,*$2,Interval::ALL_REALS);
		                                               scopes.top().add_symbol($1,*$$);  
		                                               free($1); delete $2; }
              | TK_NEW_SYMBOL dimension 
	            TK_IN const_expr                     { $$ = new P_ExprSymbol($1,*$2,*$4);
		                                               scopes.top().add_symbol($1,*$$); 
						                               free($1); delete $2; delete $4; }
              ; 

dimension     :                                      { $$=new Dim(0,0,0); }
              | '[' expr ']'                         { $$=new Dim(0,0,_2int(*$2); delete $2; }
              | '[' expr ']' '[' expr ']'            { $$=new Dim(0,_2int(*$2),_2int(*$5)); delete $2; delete $5; }
              | '[' expr ']' '[' expr ']'
 	            '[' expr ']'                         { $$=new Dim(_2int(*2),_2int(*$5),$8->_2int()); 
		                                               delete $2; delete $5; delete $8; }
	          ;

interval      : '[' expr ',' expr ']'                { $$=new Interval(_2dbl(*$2), _2dbl(*$4)); 
		                                               delete $2; delete $4; }
              ;

/**********************************************************************************************************************/
/*                                                FUNCTIONS                                                           */
/**********************************************************************************************************************/

decl_fnc_list : decl_fnc_list decl_fnc 
              | 
              ;

decl_fnc      : TK_FUNCTION                          { scopes.push(Scope(scopes.top(),TK_CONSTANT)); }
                TK_NEW_SYMBOL dimension              { scopes.top().add_func_return($3,*$4); delete $4; }
				TK_EQU TK_NEW_SYMBOL
                '(' fnc_inpt_list ')'
				fnc_code
				TK_RETURN_SYMBOL TK_EQU expr semicolon_opt
				TK_END                               { Function* f=new Function($9,$14,$3);
													   source.func.push_back(f);
													   scopes.pop();
													   scopes.top().add_func($3,f); 
                                        	       	   free($3); delete $4; free($7); }
              ;

semicolon_opt : ';' | ;

fnc_inpt_list : fnc_inpt_list ',' fnc_input          { $1->push_back(strdup($3)); free($3); $$=$1; }
              | fnc_input                            { $$=new vector<const char*>(); $$->push_back(strdup($1)); free($1);}
              ;

fnc_input     : TK_NEW_SYMBOL dimension              { scopes.top().add_func_input($1,*$2); $$=$1; delete $2; }
              ;

fnc_code      : fnc_code fnc_assign ';'              
              | fnc_assign ';'                       
              ;

fnc_assign    : TK_NEW_SYMBOL TK_EQU expr            { /* note: if this tmp symbol is not used, the expr will never be deleted */
														scopes.top().add_func_tmp_symbol($1,$3); free($1); }
              ;

/**********************************************************************************************************************/
/*                                                CONSTRAINTS                                                         */
/**********************************************************************************************************************/
decl_ctr_list : TK_CTRS
                ctr_blk_list TK_END            { $$ = $2; }
              ;
              
ctr_blk_list  : ctr_blk_list_ semi_col_opt     { $$ = $1; }
              ;

ctr_blk_list_ : ctr_blk_list_ ';' ctr_blk      { $1->add(*$3); $$ = $1; }
              | ctr_blk                        { $$ = new vector<P_NumConstraint*>(*$1); }
              ;

ctr_blk       : ctr_loop                       { $$ = $1; }
              | ctr                            { $$ = new P_OneConstraint($1); }
              ;


ctr_loop      : TK_FOR TK_NEW_SYMBOL TK_EQU
				const_expr ':' const_expr ';'  { scopes.push(scopes.top());
						       					 scopes.top().add_iterator($2); }
                ctr_blk_list TK_END            { $$ = new P_ConstraintLoop($2, $4->_2int(), $6->_2int(), *$9); 
						                         scopes.pop();
		                                         free($2); delete $4; delete $6; }
              ;

ctr           : expr TK_EQU expr               { $$ = &(*$1 =*$3); }
              | expr TK_LEQ expr               { $$ = &(*$1<=*$3); }
              | expr TK_GEQ expr               { $$ = &(*$1>=*$3); }
              | expr   '<'  expr               { $$ = &(*$1< *$3); }
              | expr   '>'  expr               { $$ = &(*$1> *$3); }
              | '(' ctr ')'                    { $$ = $2; }
              ; 
              
/**********************************************************************************************************************/
/*                                                EXPRESSIONS                                                         */
/**********************************************************************************************************************/

expr          : expr '+' expr	                     { $$ = &(*$1 + *$3); }
              | expr '*' expr	                     { $$ = &(*$1 * *$3); }
              | expr '-' expr	                     { $$ = &(*$1 - *$3); }
              | expr '/' expr	                     { $$ = &(*$1 / *$3); }
              | TK_MAX '(' expr_list ')'             { $$ = &max(*$3); delete $3; }
              | TK_MIN '(' expr_list ')'             { $$ = &min(*$3); delete $3; }
              | TK_ATAN2 '(' expr ',' expr ')'       { $$ = &atan2(*$3,*$5); }
              | '-' expr                             { $$ = &(-*$2); }
              | TK_ABS  '(' expr ')'                 { $$ = &abs  (*$3); }
              | TK_SIGN '(' expr ')'                 { $$ = &sign (*$3); }
              | expr '^' expr	                     { $$ = new P_ExprPower(*$1, *$3); }
              | TK_SQRT '(' expr ')'                 { $$ = &sqrt (*$3); }
              | TK_EXPO '(' expr ')'                 { $$ = &exp  (*$3); }
              | TK_LOG '(' expr ')'                  { $$ = &ln   (*$3); }
              | TK_COS '(' expr ')'                  { $$ = &cos  (*$3); }
              | TK_SIN '(' expr ')'                  { $$ = &sin  (*$3); }
              | TK_TAN '(' expr ')'                  { $$ = &tan  (*$3); }
              | TK_ACOS '(' expr ')'                 { $$ = &acos (*$3); }
              | TK_ASIN '(' expr ')'                 { $$ = &asin (*$3); }
              | TK_ATAN '(' expr ')'                 { $$ = &atan (*$3); }
              | TK_COSH '(' expr ')'                 { $$ = &cosh (*$3); }
              | TK_SINH '(' expr ')'                 { $$ = &sinh (*$3); }
              | TK_TANH '(' expr ')'                 { $$ = &tanh (*$3); }
              | TK_ACOSH '(' expr ')'                { $$ = &acosh(*$3); }
              | TK_ASINH '(' expr ')'                { $$ = &asinh(*$3); }
              | TK_ATANH '(' expr ')'                { $$ = &atanh(*$3); }
              | '+' expr                             { $$ = $2; }
              | '(' expr ')'		                 { $$ = $2; }
              | TK_SYMBOL                            { $$ = scope.get_symbol($1); free($1); }
              | TK_TMP_SYMBOL                        { $$ = scope.get_expr($1); free($1); }
              | TK_CONSTANT                          { $$ = new ExprConstant(scope.get_cst($1)); free($1); }
              | TK_FUNC '(' expr_list ')'            { $$ = scope.get_func($1)(*$3); free($1); delete $3; }
              | '[' expr_list ']'                    { $$ = new ExprVector(*$2); delete $2;
              | expr '[' expr ']'                    { $$ = new P_ExprIndex(*$1,*$3); }
              | TK_FLOAT                             { $$ = new ExprConstant($1); }
              | TK_INTEGER                           { $$ = new ExprConstant($1); }
              | interval                             { $$ = new ExprConstant($1); delete $1; }
	      	  ;
	      
expr_list     : expr_list ',' expr                   { $1->push_back($3); $$=$1; }
              | expr                                 { $$ = new vector<const ExprNode*>; 
              										   $$->push_back($1); }
              ;