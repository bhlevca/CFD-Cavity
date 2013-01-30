/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef TOKENIZER_INCLUDED
#define TOKENIZER_INCLUDED
#include <stdio.h>
#define MAXLINELENGTH 32768
#define TOKEN_EOF -1

typedef struct {
    const char* token;
    int         tokenid;
} Token;


//! Class Tokenizer contains some methods important for the NAV-file parser 
class Tokenizer {
public:
     //! read a value of type char* q=new char[..]
    char*       GetStringBB(char b);
    void        SetEndTokenChars(char* pc);
    //! read a new line
    void        NewLine();
    //! read a value of type double
    double      GetDouble();
    //! read a value of type long
    long        GetLong();
    //! read a value of type char[..] 
    char*       GetString();
    //! read next Token 
    int         GetToken(const Token* token_arr,int jumpback=0);
    int         ShowToken(const Token* token_arr);
               ~Tokenizer();
                Tokenizer(FILE* file);
                Tokenizer(const char* filename);
    void        Error(const char* err);

protected:
    int         GetNextItem();
    FILE*       f;
    int         linenumber;
    char        actual_line[MAXLINELENGTH];
    char*       actual_ptr;
    int         toka,toka_token;
    char*       endtokenchars;
};

#endif
