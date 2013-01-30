/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#include "tokenizer.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

char* etchar="";

Tokenizer::Tokenizer(FILE* file) {
    f=file;assert(f);
    linenumber=0;actual_line[0]=0;actual_ptr=actual_line;
    toka=0;endtokenchars=etchar;
}

Tokenizer::Tokenizer(const char* filename) {
    f=fopen(filename,"r");
    if(!f) {
        fprintf(stderr,"Tokenizer: Cannot open file %s.\n",filename);
        exit(1);
    }
    linenumber=0;actual_line[0]=0;actual_ptr=actual_line;
    toka=0;endtokenchars=etchar;
}

Tokenizer::~Tokenizer() {
    if(f) fclose(f);
}

void Tokenizer::Error(const char* err) {
    fprintf(stderr,"Error in line %d (c%d): %s:\n\n%s",linenumber,(int)(actual_ptr-actual_line)+1,err,actual_line);
    for(int i=0;i<(int)(actual_ptr-actual_line)-1;i++) fprintf(stderr," ");
    fprintf(stderr,"^^^\n\n");
    exit(1);
}

int Tokenizer::GetToken(const Token * token_arr,int jumpback) {
    if((*actual_ptr==0) && feof(f)) return TOKEN_EOF;
    if(!GetNextItem()) return TOKEN_EOF;
    while(token_arr->tokenid!=0) {
        int sl=strlen(token_arr->token);
        if(strncmp(token_arr->token,actual_ptr,sl)==0 && 
            (isspace(*(actual_ptr+sl)) || strchr(endtokenchars,*(actual_ptr+sl))!=0 || strlen(token_arr->token)==1)) {
            if(!jumpback) actual_ptr+=sl;
            return token_arr->tokenid;
        }
        token_arr++;
    }
    Error("Syntax Error");
    return TOKEN_EOF;
}

int Tokenizer::ShowToken(const Token* token_arr) {
    toka_token=GetToken(token_arr,1);   
    return toka_token;
}

int Tokenizer::GetNextItem() {
    do {
        if((*actual_ptr==0) && feof(f)) return 0;
        if(*actual_ptr==0) if(!fgets(actual_line,MAXLINELENGTH,f)) return 0; else {
            linenumber++;
            actual_ptr=actual_line;
        }
        while(isspace(*actual_ptr)) actual_ptr++;
    } while(*actual_ptr==0);
    return 1;
}

char* Tokenizer::GetString() {
    if((*actual_ptr==0) && feof(f)) return 0;
    if(!GetNextItem()) Error("Premature end of file.");
    char* p=actual_ptr,*q;
    while(!isspace(*p) && *p!=0) p++;
    int count=(int)(p-actual_ptr);
    q=new char[count+1];
    for(int i=0;i<count;i++) q[i]=actual_ptr[i];
    q[count]=0;
    actual_ptr=p;
    toka=0;
    return q;
}

double Tokenizer::GetDouble() {
    if((*actual_ptr==0) && feof(f)) Error("End of file, DOUBLE expected");
    if(!GetNextItem()) Error("Premature end of file, DOUBLE expected");
    char* p;
    double d=strtod(actual_ptr,&p);
    if(p==actual_ptr) Error("Cannot parse DOUBLE");
    actual_ptr=p;
    toka=0;
    return d;
}

long Tokenizer::GetLong() {
    if((*actual_ptr==0) && feof(f)) Error("End of file, DOUBLE expected");
    if(!GetNextItem()) Error("Premature end of file, DOUBLE expected");
    char* p;
    long d=strtol(actual_ptr,&p,0);
    if(p==actual_ptr) Error("Cannot parse LONG");
    actual_ptr=p;
    toka=0;
    return d;
}

void Tokenizer::NewLine() {
    *actual_ptr=0;
//  GetNextItem();
}

void Tokenizer::SetEndTokenChars(char * pc) {
    endtokenchars=pc;
}

char* Tokenizer::GetStringBB(char b) {
    while((*actual_ptr!=b || isspace(*actual_ptr)) && *actual_ptr!=0) actual_ptr++;
    if(*actual_ptr==0) return 0;
    actual_ptr++;
    char* p=strchr(actual_ptr,b);
    if(!p) return 0;
    char* q=new char[(int)(p-actual_ptr)+1];
    strncpy(q,actual_ptr,(int)(p-actual_ptr));
    q[(int)(p-actual_ptr)]=0;
    actual_ptr=p+1;
    toka=0;
    return q;
}
