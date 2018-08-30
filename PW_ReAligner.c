#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#include <stdint.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#define Max_Bandwidth 2000lu
#define Max_MA_Breadth 700000lu
#define Max_Seq_Length 35000lu
#define Max_Seq_Anzahl 18000lu
#define Max_Long 18446744073709551615lu

//#define PRINT

//#define CURRENTFUNCTION

/*
This tool realignes sequences to a multiple sequence alignment until it converges. 
It optimizes the sum of all induced pairwise-alignment-scores.
The input is a msa with equally long rows of 'a','c','g','t','-',' ' ended by '\n'.
*/

unsigned long Matrix[Max_Seq_Length][Max_Bandwidth]; 
int Way[Max_Seq_Length];  
char Seq_Bases[Max_Seq_Length];
int bandwidth;
int bandwidthhalf;
int length;
int Gaps[120];
char Chars[6];
int Lengths[Max_Seq_Anzahl];

//The data structure of a column: These are combined into a linked list which represents the msa
typedef struct Column
{
	char *Bases;
	struct Column* Previous;
	struct Column* Next;
	unsigned long *w_con;
}Column;

//The columns are accessed via a weighted consensus array that is used for the alignment
Column *W_Con2Columns[Max_MA_Breadth];
Column *Reservoir[Max_MA_Breadth];
int reserve;


void *Guarded_Malloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL)
    { fprintf(stderr,"\nError in PW_ReAligner\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}

void *Guarded_Realloc(void *p, size_t size)
{ p = realloc(p,size);
  if (p == NULL)
    { fprintf(stderr,"\nError in PW_ReAligner:\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}  

static inline unsigned long u_intmin(unsigned long x, unsigned long y){if(x<y){return x;} return y;}

static inline unsigned long u_intmax(unsigned long x, unsigned long y){if(x>y){return x;} return y;}

static inline int intmin(int x, int y){if(x<y){return x;} return y;}

static inline int intmax(int x, int y){if(x>y){return x;} return y;}

/******************************************MMA_Einlesen*****************************************/
int Tiefe;
int Breite;
unsigned long** w_con_array;
Column* Null_Column;
Column* Last_Column;

//reading in the msa
void MMA_Einlesen(char *MApath)
{		
  #ifdef CURRENTFUNCTION
  printf("In MMA_Einlesen\n");
  fflush(stdout);
  #endif

  char buffer[Max_MA_Breadth];
  register char *s;
  int i;
  w_con_array=Guarded_Malloc(sizeof(unsigned long*)*Max_MA_Breadth);  //Guarded_Malloc(sizeof(char*)*Max_MA_Breadth);
  Column* Present_Column;
  Column* PreviousColumn;
  Null_Column=Guarded_Malloc(sizeof(Column));
  Null_Column->Bases=Guarded_Malloc(sizeof(char)*Max_Seq_Anzahl);
  Null_Column->w_con=Guarded_Malloc(sizeof(unsigned long)*6);
  *(Null_Column->w_con)=0;
  *(Null_Column->w_con+1)=0;
  *(Null_Column->w_con+2)=0;
  *(Null_Column->w_con+3)=0;
  *(Null_Column->w_con+4)=0;
  *(Null_Column->w_con+5)=0;

  int Reihe=-1;
  //while ((s=fgets(buffer, Max_MA_Breadth-2, stdin)) != NULL)
  FILE * File;
  size_t len=Max_MA_Breadth-2;
  File=fopen(MApath,"r");
  if(File==NULL){printf("MA is missing.\n"); exit(1);}
  while ((s = fgets(buffer, len, File)) != NULL)  
  {

    Reihe++;

    #ifdef PRINT
    printf("%d\n",Reihe);
    fflush(stdout);
    #endif

    Breite = strlen(buffer);

    if (buffer[Breite-1] != '\n'){printf("Was für einen Scheiß liest Du da ein?\n");exit(1);}

    Breite--;  //length of the row included '\n'

    if (Reihe==0)  //allocating memory for the msa
    {
      Present_Column=Null_Column;
      for(i=1;i<Breite;i++)
      {
      	Present_Column->Next=Guarded_Malloc(sizeof(Column));
      	PreviousColumn=Present_Column;
      	Present_Column=Present_Column->Next;
      	Present_Column->Bases=Guarded_Malloc(sizeof(char)*Max_Seq_Anzahl);
      	Present_Column->w_con=Guarded_Malloc(sizeof(unsigned long)*6);
		    *(Present_Column->w_con)=0;
	      *(Present_Column->w_con+1)=0;
    	  *(Present_Column->w_con+2)=0;
		    *(Present_Column->w_con+3)=0;
		    *(Present_Column->w_con+4)=0;
		    *(Present_Column->w_con+5)=0;      	
      	Present_Column->Previous=PreviousColumn;
      }
      Last_Column=Present_Column;
    }

    Lengths[Reihe]=0;
    //Entering the bases.
    Present_Column=Null_Column;
    for(i=0;i<Breite;i++)
    {

      if(buffer[i]=='a' || buffer[i]=='A')
      {
        Lengths[Reihe]++;
      	*(Present_Column->Bases+Reihe)=0;

	      *(Present_Column->w_con+1)+=1;
    	  *(Present_Column->w_con+2)+=1;
    		*(Present_Column->w_con+3)+=1;
    		*(Present_Column->w_con+4)+=1;
    		*(Present_Column->w_con+5)+=1;            	
      }
      else if(buffer[i]=='c' || buffer[i]=='C')
      {
        Lengths[Reihe]++;
        *(Present_Column->Bases+Reihe)=1;
		    *(Present_Column->w_con)+=1;

    	  *(Present_Column->w_con+2)+=1;
		    *(Present_Column->w_con+3)+=1;
		    *(Present_Column->w_con+4)+=1;
		    *(Present_Column->w_con+5)+=1;         
      }
      else if(buffer[i]=='g' || buffer[i]=='G')
      {
        Lengths[Reihe]++;
      	*(Present_Column->Bases+Reihe)=2;
		    *(Present_Column->w_con)+=1;
	      *(Present_Column->w_con+1)+=1;

		    *(Present_Column->w_con+3)+=1;
		    *(Present_Column->w_con+4)+=1;
		    *(Present_Column->w_con+5)+=1;       	
      }
      else if(buffer[i]=='t' || buffer[i]=='T')
      {
        Lengths[Reihe]++;
      	*(Present_Column->Bases+Reihe)=3;
		    *(Present_Column->w_con)+=1;
	      *(Present_Column->w_con+1)+=1;
    	  *(Present_Column->w_con+2)+=1;

		    *(Present_Column->w_con+4)+=1;
		    *(Present_Column->w_con+5)+=1;       	
      }
      else if(buffer[i]=='-' || buffer[i]=='_')
      {
      	*(Present_Column->Bases+Reihe)=4;
		    *(Present_Column->w_con)+=1;
	      *(Present_Column->w_con+1)+=1;
    	  *(Present_Column->w_con+2)+=1;
		    *(Present_Column->w_con+3)+=1;

		    *(Present_Column->w_con+5)+=1;       	
      }
      else if(buffer[i]==' ')
      {
      	*(Present_Column->Bases+Reihe)=5;
      }
      Present_Column=Present_Column->Next;
    }
  }
  Tiefe=Reihe+1;

  //Allocating a small reservoir of columns: 
  reserve=0;
  for(i=0;i<Breite/10;i++)
  {
    reserve++;
    Reservoir[i]=Guarded_Malloc(sizeof(Column));
    Reservoir[i]->Bases=Guarded_Malloc(sizeof(char)*Max_Seq_Anzahl);
    Reservoir[i]->w_con=Guarded_Malloc(sizeof(unsigned long)*6);    
  }
  #ifdef CURRENTFUNCTION
  printf("Out MMA_Einlesen\n");
  fflush(stdout);
  #endif
}

static inline unsigned long Score(int index,char Base)
{
  return *(*(w_con_array+index)+Base);
}

//Reading from the alignment matrix
static inline unsigned long MatrixOut(int x, int y)
{
  #ifdef CURRENTFUNCTION
  printf("In MatrixOut\n");
  fflush(stdout);
  #endif

  if(x==-1)
  {
    #ifdef CURRENTFUNCTION
    printf("Out MatrixOutx-1 %d\n",y);
    fflush(stdout);
    #endif    
    return 0;
  }  //Before the start of a sequence ' ' can be added.

  if(y==-1)
  {
  #ifdef CURRENTFUNCTION
  printf("Out MatrixOuty-1\n");
  fflush(stdout);
  #endif
  return Max_Long/2;
  }  //Max_MA_Breadth*Max_Seq_Anzahl;}    //No aligning beyond the msa.

  int anf=intmax(0,Way[x]-bandwidthhalf);

  if(y-anf<0)  //Accessing a non calculated matrix entry
  {
    #ifdef CURRENTFUNCTION
    printf("Out MatrixOut<0\n");
    fflush(stdout);
    #endif
    return Max_Long/2;
  }

  if(y-anf>bandwidth-1)  //This handles big gaps. We possibly have to jump through a non-calculated part of the matrix
  {
    if(x==length-1)return Matrix[x][bandwidth-1];   //Adding ' '
    unsigned long score=Matrix[x][bandwidth-1];
    while(y-anf>bandwidth-1){score+=Score(y,4);y--;}
    #ifdef CURRENTFUNCTION
    printf("Out MatrixOutrec\n");
    fflush(stdout);
    #endif
    return score;
  }

  #ifdef CURRENTFUNCTION
  printf("Out MatrixOut\n");
  fflush(stdout);
  #endif

  return Matrix[x][y-anf];
}

// Writing into the alignment matrix
static inline void MatrixIn(int x, int y, unsigned long input)
{
  #ifdef CURRENTFUNCTION
  printf("In MatrixIn\n");
  fflush(stdout);
  #endif

  if(y==-1){y=Max_Bandwidth-1;}
  if(x==-1){x=Max_Seq_Length-1;}

  int anf=intmax(0,Way[x]-bandwidthhalf);
  Matrix[x][y-anf]=input;	

  #ifdef CURRENTFUNCTION
  printf("Out MatrixIn\n");
  fflush(stdout);
  #endif
}



//'-' is turned into ' '
int EntAlGapper2()
{

  printf("Enter Algapper\n");
  fflush(stdout);

  int i,k,count,count2;
  Column* Present_Column;
  Column* Speicher;

  Present_Column=Null_Column;

  count=0;
  count2=0;

  int basevorhanden=0;

  for(i=0;i<Breite;i++)
  {
    printf("%d - %d - %d\n",i,Breite,count2 );
    fflush(stdout);

    basevorhanden=0;
    for(k=0;k<Tiefe;k++)
    {

      // At the start there can be no alignment gap. 
      if(i==0 && *(Present_Column->Bases+k)==4)
      {
        *(Present_Column->Bases+k)=5;
        
        *(Present_Column->w_con+0)-=1;
        *(Present_Column->w_con+1)-=1;
        *(Present_Column->w_con+2)-=1;
        *(Present_Column->w_con+3)-=1;
        *(Present_Column->w_con+5)-=1;
        count++;        
      }

      else if(i>0 && *(Present_Column->Bases+k)==4 && *((Present_Column->Previous)->Bases+k)==5)
      {
        *(Present_Column->Bases+k)=5;
        
        *(Present_Column->w_con+0)-=1;
        *(Present_Column->w_con+1)-=1;
        *(Present_Column->w_con+2)-=1;
        *(Present_Column->w_con+3)-=1;
        *(Present_Column->w_con+5)-=1;
        count++;
      }

      if(*(Present_Column->Bases+k)<4)basevorhanden=1;
    }

    Speicher=Present_Column->Next;

    if(!basevorhanden)  //Columns without base are deleted
    {
      printf("Inderlöscherei\n");
      fflush(stdout);

      count2++;
      //Extreme cases first: 
      if(Present_Column==Null_Column)
      {
        printf("if\n");
        fflush(stdout);
        Null_Column=Present_Column->Next;
        free(Present_Column->Bases);
        free(Present_Column->w_con);
        free(Present_Column);
        Breite--;
        i--;
      }
      else if(Present_Column==Last_Column)
      {
        printf("else if\n");
        fflush(stdout);        
        Last_Column=Present_Column->Previous;
        free(Present_Column->Bases);
        free(Present_Column->w_con);
        free(Present_Column);
        Breite--;  
        i--;      
      }
      else
      {
        printf("else\n");
        fflush(stdout);        
        (Present_Column->Previous)->Next=Present_Column->Next;
        printf("else1\n");
        fflush(stdout);          
        (Present_Column->Next)->Previous=Present_Column->Previous;
        printf("else2\n");
        fflush(stdout);  
        
        free(Present_Column->Bases);

        printf("else3\n");
        fflush(stdout);          
        free(Present_Column->w_con);
        printf("else4\n");
        fflush(stdout);          
        free(Present_Column);
        printf("else5\n");
        fflush(stdout);          
        Breite--;     
        i--;      
      }

      printf("deleted\n");
      fflush(stdout);
    }

    Present_Column=Speicher;
  }  

  printf("Leave Algapper\n");
  fflush(stdout);
  return count;



/*  #ifdef PRINT
  printf("%d AlGaps zu CovGaps umgewandelt und %d Columns gelöscht.\n",count,count2);  
  #endif
  return count;*/
}


//Turning '-' into ' '
int EntAlGapper()
{
  #ifdef CURRENTFUNCTION
  printf("In EntAlGapper\n");
  fflush(stdout);
  #endif

  int i,k,count,count2;
  Column* Present_Column;
  Column* Speicher;

  Present_Column=Null_Column;
  count=0;
  count2=0;
  int basevorhanden=0;

  for(i=0;i<Breite;i++)
  {

    basevorhanden=0;
    for(k=0;k<Tiefe;k++)
    {

      // Before the sequence starts there can be no alignment gap
      if(i==0 && *(Present_Column->Bases+k)==4)
      {
        *(Present_Column->Bases+k)=5;
        
        *(Present_Column->w_con+0)-=1;
        *(Present_Column->w_con+1)-=1;
        *(Present_Column->w_con+2)-=1;
        *(Present_Column->w_con+3)-=1;
        *(Present_Column->w_con+5)-=1;
        count++;        
      }

      else if(i>0 && *(Present_Column->Bases+k)==4 && *(W_Con2Columns[i-1]->Bases+k)==5)
      {
        *(Present_Column->Bases+k)=5;
        
        *(Present_Column->w_con+0)-=1;
        *(Present_Column->w_con+1)-=1;
        *(Present_Column->w_con+2)-=1;
        *(Present_Column->w_con+3)-=1;
        *(Present_Column->w_con+5)-=1;
        count++;
      }

      if(*(Present_Column->Bases+k)<4)basevorhanden=1;
    }

    Speicher=Present_Column->Next;

    //To be able to run backwards through the columns
    W_Con2Columns[i]=Present_Column;

    if(!basevorhanden)  //A column without bases is deleted
    {
      count2++;
      //Extreme cases
      if(Present_Column==Null_Column)
      {
        Null_Column=Present_Column->Next;

        Reservoir[reserve]=Present_Column;
        reserve++;

        Breite--;
        i--;
      }
      else if(Present_Column==Last_Column)
      {
        Last_Column=Present_Column->Previous;

        Reservoir[reserve]=Present_Column;
        reserve++;
        Breite--;  
        i--;      
      }
      else
      {
        (Present_Column->Previous)->Next=Present_Column->Next;
        (Present_Column->Next)->Previous=Present_Column->Previous;

        Reservoir[reserve]=Present_Column;
        reserve++;

        Breite--;     
        i--;      
      }
    }

    Present_Column=Speicher;
  }  


  // Running backwards through the columns 
  for(i=Breite-1;i>=0;i--)
  {

    Present_Column=W_Con2Columns[i];

    basevorhanden=0;
    for(k=0;k<Tiefe;k++)
    {

      //No alignment gaps at the end of a sequence
      if(i==Breite-1 && *(Present_Column->Bases+k)==4)
      {
        *(Present_Column->Bases+k)=5;
        
        *(Present_Column->w_con+0)-=1;
        *(Present_Column->w_con+1)-=1;
        *(Present_Column->w_con+2)-=1;
        *(Present_Column->w_con+3)-=1;
        *(Present_Column->w_con+5)-=1;
        count++;        
      }

      else if(i<Breite-1 && *(Present_Column->Bases+k)==4 && *(W_Con2Columns[i+1]->Bases+k)==5)
      {
        *(Present_Column->Bases+k)=5;
        
        *(Present_Column->w_con+0)-=1;
        *(Present_Column->w_con+1)-=1;
        *(Present_Column->w_con+2)-=1;
        *(Present_Column->w_con+3)-=1;
        *(Present_Column->w_con+5)-=1;
        count++;
      }  

      if(*(Present_Column->Bases+k)<4)basevorhanden=1;
    }

  
    if(!basevorhanden)  //Columns without bases are deleted
    {
      count2++;
      if(Present_Column==Null_Column)
      {
        Null_Column=Present_Column->Next;
/*        free(Present_Column->Bases);
        free(Present_Column->w_con);
        free(Present_Column);*/
        //Ins Reservoir statt freeen.
        Reservoir[reserve]=Present_Column;
        reserve++;

        Breite--;
      }
      else if(Present_Column==Last_Column)
      {
        Last_Column=Present_Column->Previous;
/*        free(Present_Column->Bases);
        free(Present_Column->w_con);
        free(Present_Column);*/
        Reservoir[reserve]=Present_Column;
        reserve++;

        Breite--;      
      }
      else
      {
        (Present_Column->Previous)->Next=Present_Column->Next;
        (Present_Column->Next)->Previous=Present_Column->Previous;

/*        free(Present_Column->Bases);
        free(Present_Column->w_con);
        free(Present_Column);*/
        Reservoir[reserve]=Present_Column;
        reserve++;

        Breite--;         
      }
    }
  } 

  #ifdef CURRENTFUNCTION
  printf("Out EntAlGapper\n");
  fflush(stdout);
  #endif

  #ifdef PRINT
  printf("%d algaps turned into covgaps %d columns deleted.\n",count,count2);  
  #endif
  return count;
}

int TheWay(int k)
{

  int i,count,lag;
  int gapcount=0;
  Column* Present_Column;
  Present_Column=Null_Column;
  count=0;
  lag=0;
  for(i=0;i<Breite;i++)
  {
 
    if(i<Breite-1 && *(Present_Column->Bases+k)==4 && *((Present_Column->Next)->Bases+k)==5)
    {
      #ifdef PRINT
      printf("Algap next to CovGap! %d-%d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",count,i);
      #endif
      //exit(1);
    }

    if(i<Breite-1 && *(Present_Column->Bases+k)<4 && *((Present_Column->Next)->Bases+k)==5)
    {
      Gaps[gapcount]=count; //The last gap before a coverage gap
      gapcount++;
    }

    if(*(Present_Column->Bases+k)<4)
    {  
      if(count>Max_Seq_Length)
      {
        printf("Max_Seq_Length to small.%d\n",count);
        fflush(stdout);
        exit(1);
      }
      Way[count]=i-lag;
/*      if(*(Present_Column->w_con+4)==1)
      {
        lag++;
      } */ 
      Seq_Bases[count]=*(Present_Column->Bases+k);  
      count++;
    }

    Present_Column=Present_Column->Next;
  }

  if(Gaps[gapcount-1]!=count-1)
  {
    Gaps[gapcount]=count-1;   //gap at the end ...
    Gaps[gapcount+1]=Max_Seq_Length; //the end
  }
  else
  {
    Gaps[gapcount]=Max_Seq_Length;
  }
  if(count!=Lengths[k])printf("TheWay miscounted %d %d\n", count,Lengths[k]);

  return count;  
}
void W_Con()
{
  #ifdef CURRENTFUNCTION
  printf("In W_Con\n");
  fflush(stdout);
  #endif  

  int i;
  Column* Present_Column;
  Column* Speicher;
  Present_Column=Null_Column;
  i=0;
  while(i<Breite)
  {
  	*(w_con_array+i)=Present_Column->w_con;

    Speicher=Present_Column->Next;

    if(*(Present_Column->w_con+4)==0)  //Deleting columns without base
    {
      if(Present_Column==Null_Column)
      {
        Null_Column=Present_Column->Next;
        Reservoir[reserve]=Present_Column;
        reserve++;
        Breite--;
      }
      else if(Present_Column==Last_Column)
      {
        Last_Column=Present_Column->Previous;
        Reservoir[reserve]=Present_Column;
        reserve++;
        Breite--;        
      }
      else
      {
        (Present_Column->Previous)->Next=Present_Column->Next;
        (Present_Column->Next)->Previous=Present_Column->Previous;

        Reservoir[reserve]=Present_Column;
        reserve++;
        Breite--;           
      }
    }

    else  //if the column has been deleted we don't progress down the w_con_array
    {
      i++;
    }

    Present_Column=Speicher;
  }

  #ifdef CURRENTFUNCTION
  printf("Out W_Con\n");
  fflush(stdout);
  #endif   
}

int W_Con_Checker()
{
  #ifdef PRINT
  printf("W_Con_Checker started.\n");
  fflush(stdout);
  #endif
  int i,Reihe,count;
  Column* Present_Column;
  Present_Column=Null_Column;
  *w_con_array=Present_Column->w_con;
  unsigned long *test_w_con;
  test_w_con=Guarded_Malloc(sizeof(unsigned long)*6);
  count=0;

  for(i=0;i<Breite-1;i++)
  {
    *(test_w_con)=0;
    *(test_w_con+1)=0;
    *(test_w_con+2)=0;
    *(test_w_con+3)=0;
    *(test_w_con+4)=0;
    *(test_w_con+5)=0;

    for(Reihe=0;Reihe<Tiefe;Reihe++)
    {     
      if(*(Present_Column->Bases+Reihe)==0)
      {

        *(test_w_con+1)+=1;
        *(test_w_con+2)+=1;
        *(test_w_con+3)+=1;
        *(test_w_con+4)+=1;
        *(test_w_con+5)+=1;              
      }
      else if(*(Present_Column->Bases+Reihe)==1)
      {
        *(test_w_con)+=1;

        *(test_w_con+2)+=1;
        *(test_w_con+3)+=1;
        *(test_w_con+4)+=1;
        *(test_w_con+5)+=1;         
      }
      else if(*(Present_Column->Bases+Reihe)==2)
      {
        *(test_w_con)+=1;
        *(test_w_con+1)+=1;

        *(test_w_con+3)+=1;
        *(test_w_con+4)+=1;
        *(test_w_con+5)+=1;        
      }
      else if(*(Present_Column->Bases+Reihe)==3)
      {
        *(test_w_con)+=1;
        *(test_w_con+1)+=1;
        *(test_w_con+2)+=1;

        *(test_w_con+4)+=1;
        *(test_w_con+5)+=1;        
      }
      else if(*(Present_Column->Bases+Reihe)==4)
      {
        *(test_w_con)+=1;
        *(test_w_con+1)+=1;
        *(test_w_con+2)+=1;
        *(test_w_con+3)+=1;

        *(test_w_con+5)+=1;        
      }
    }

    //Checking the weighted consensus:
    int ok;
    ok=1;
    if(*(test_w_con) != *(*(w_con_array+i)))
    {count++;printf("w_con_gau0, column %d!!!!%lu - %lu\n",i,*(test_w_con),*(*(w_con_array+i)));ok=0;}
    if(*(test_w_con+1) != *(*(w_con_array+i)+1))
    {count++;printf("w_con_gau1, column %d!!!!%lu - %lu!!!!\n",i,*(test_w_con+1),*(*(w_con_array+i)+1));ok=0;}
    if(*(test_w_con+2) != *(*(w_con_array+i)+2))
    {count++;printf("w_con_gau2, column %d!!!!%lu - %lu!!!!\n",i,*(test_w_con+2),*(*(w_con_array+i)+2));ok=0;}
    if(*(test_w_con+3) != *(*(w_con_array+i)+3))
    {count++;printf("w_con_gau3, column %d!!!!%lu - %lu!!!!\n",i,*(test_w_con+3),*(*(w_con_array+i)+3));ok=0;}
    if(*(test_w_con+4) != *(*(w_con_array+i)+4))
    {count++;printf("w_con_gau4, column %d!!!!%lu - %lu!!!!\n",i,*(test_w_con+4),*(*(w_con_array+i)+4));ok=0;}
    if(*(test_w_con+5) != *(*(w_con_array+i)+5))
    {count++;printf("w_con_gau5, column %d!!!!%lu - %lu!!!!\n",i,*(test_w_con+5),*(*(w_con_array+i)+5));ok=0;}
    if(!ok)exit(1);

    //Next column:
    Present_Column=Present_Column->Next;
  }
  free(test_w_con);
  return count;
}




unsigned long ReihenScore(int k)
{
  #ifdef CURRENTFUNCTION
  printf("In ReihenScore\n");
  fflush(stdout);
  #endif 

  int i;
  unsigned long score=0;
  //length=TheWay(k);
  W_Con();
  Column* Present_Column;
  Present_Column=Null_Column;
  char Base;
  for(i=0;i<Breite;i++)
  {
    Base=*(Present_Column->Bases+k);
    if(Base!=5){score+=Score(i,Base);}
    Present_Column=Present_Column->Next;
  }
  
  //printf("%d,%03d,%03d\n",score/1000000,(score/1000)%1000,score%1000 );
  #ifdef CURRENTFUNCTION
  printf("Out ReihenScore\n");
  fflush(stdout);
  #endif 

  return score;
}

unsigned long OverallScore()
{
  int k;
  unsigned long score=0;

  W_Con();
  for(k=0;k<Tiefe;k++)
  {
    score+=ReihenScore(k);
  }
  return score;
}

//buggy version? 
unsigned long OverallScore2()
{
  int i,j,k;
  unsigned long score=0;
  int bases[5];
  W_Con();
  for(i=0;i<Breite;i++)
  {
    k=0;
    for(j=0;j<5;j++)k+=Score(i,j);
    k/=4; //The Coverage
    for(j=0;j<5;j++)bases[j]=k-Score(i,j); //how many bases
    for(j=0;j<5;j++)
    {
      for(k=0;k<5;k++)
      {
        if(k!=j)score+=bases[j]*bases[k];
      }
    }
  }
  return score;
}

unsigned long BestMille;
unsigned long BestUno;
int OverallScorePrint()
{
  //int i,j,k;
  int k;
  unsigned long scoreMille=0;
  unsigned long scoreUno=0;
  //int bases[5];
  W_Con();

  for(k=0;k<Tiefe;k++)
  {
    scoreUno+=ReihenScore(k);
    while(scoreUno>1000000)
    {
      scoreUno-=1000000;
      scoreMille+=1;
    }
  }


  int returnint=0;
  if(scoreMille<BestMille || (scoreMille==BestMille && scoreUno<BestUno))
  {
    returnint=1;
    BestMille=scoreMille;
    BestUno=scoreUno;
  }

  printf("OverallScore: %lu%06lu\n",scoreMille,scoreUno);
  return returnint;
}  

void OverallScorePrintF(char *path, int alignmentcount)
{
  FILE *datei;
  datei=fopen(path,"a");

  if(NULL == datei)
  {
    printf("DateiVerbratei!\n");
    exit(1);
  }

  //int i,j,k;
  int k;
  unsigned long scoreMille=0;
  unsigned long scoreUno=0;
  //int bases[5];
  W_Con();
/*  for(i=0;i<Breite;i++)
  {
    k=0;
    for(j=0;j<5;j++)k+=Score(i,j);
    k/=4; //Die Coverage
    for(j=0;j<5;j++)bases[j]=k-Score(i,j); //wie viele Basen
    for(j=0;j<5;j++)
    {
      for(k=0;k<5;k++)
      {
        if(k!=j)scoreUno+=bases[j]*bases[k];
      }
    }*/
  for(k=0;k<Tiefe;k++)
  {
    scoreUno+=ReihenScore(k);
    while(scoreUno>1000000)
    {
      scoreUno-=1000000;
      scoreMille+=1;
    }
  }

  fprintf(datei,"%d %lu%06lu %d\n", alignmentcount,scoreMille,scoreUno,Breite);
  fclose(datei); 
}  

// Experimental code
double OverallQuality()
{
  double quality;
  int i,j,k;
  long int score=0;
  long int cov2=0;
  int bases[5];
  W_Con();
  for(i=0;i<Breite;i++)
  {
    k=0;
    for(j=0;j<5;j++)k+=Score(i,j);
    k/=4; //Die Coverage
    cov2+=Score(i,4)*k;  //all comparisons
    for(j=0;j<5;j++)bases[j]=k-Score(i,j); //bases

    for(j=0;j<4;j++)  //bases
    {
      for(k=0;k<5;k++)
      {
        if(k!=j)score+=bases[j]*bases[k];  //errors
      }
    }
  }
  int overalllength=0;
  for(i=0;i<Tiefe;i++)overalllength+=TheWay(i);
  overalllength/=Tiefe;

  quality=(double)score;

  //quality/=(double)(overalllength*Tiefe*Tiefe);
  quality/=(double)cov2;

  return quality;
}

// Orthogonal optimizer, probably doesn't add much
void Double_Column_Optimizer(Column* Present_Column)
{
  Column* Next_Column=Present_Column->Next;
  int i,ii,iii,iiii;
  int j,jj;
  int gapbase[4];
  gapbase[0]=0;
  gapbase[1]=0;
  gapbase[2]=0;
  gapbase[3]=0;
  unsigned long w_con1[6];
  for(i=0;i<6;i++)w_con1[i]=0;
  int w_con2[6];
  for(i=0;i<6;i++)w_con2[i]=0;

  for(i=0;i<Tiefe;i++)
  {
    if( (*(Present_Column->Bases+i)==0 && *(Next_Column->Bases+i)==4) || (*(Present_Column->Bases+i)==4 && *(Next_Column->Bases+i)==0) )gapbase[0]++;
    else if( (*(Present_Column->Bases+i)==1 && *(Next_Column->Bases+i)==4) || (*(Present_Column->Bases+i)==4 && *(Next_Column->Bases+i)==1) )gapbase[1]++;
    else if( (*(Present_Column->Bases+i)==2 && *(Next_Column->Bases+i)==4) || (*(Present_Column->Bases+i)==4 && *(Next_Column->Bases+i)==2) )gapbase[2]++;
    else if( (*(Present_Column->Bases+i)==3 && *(Next_Column->Bases+i)==4) || (*(Present_Column->Bases+i)==4 && *(Next_Column->Bases+i)==3) )gapbase[3]++;
    else 
    {
      for(ii=0;ii<6;ii++)
      {
        if(*(Present_Column->Bases+i)!=ii && *(Present_Column->Bases+i)!=5){w_con1[ii]++;}
        if(*(Next_Column->Bases+i)!=ii && *(Next_Column->Bases+i)!=5){w_con2[ii]++;}
      }
    }
  }

  int min=INT_MAX;

  int is[4];
  int is2[4];
  unsigned long score;
  for(i=0;i<2;i++)
  {
    for(ii=0;ii<2;ii++)
    {
      for(iii=0;iii<2;iii++)
      {
        for(iiii=0;iiii<2;iiii++)
        {
          score=0;

          //score with the rest
          if(i){score+=(w_con1[0]+w_con2[4])*gapbase[0];}
          else{score+=(w_con2[0]+w_con1[4])*gapbase[0];}

          if(ii){score+=(w_con1[1]+w_con2[4])*gapbase[1];}
          else{score+=(w_con2[1]+w_con1[4])*gapbase[1];}

          if(iii){score+=(w_con1[2]+w_con2[4])*gapbase[2];}
          else{score+=(w_con2[2]+w_con1[4])*gapbase[2];}

          if(iiii){score+=(w_con1[3]+w_con2[4])*gapbase[3];}
          else{score+=(w_con2[3]+w_con1[4])*gapbase[3];}   

          //score among each other: 
          is2[0]=i;
          is2[1]=ii;
          is2[2]=iii;
          is2[3]=iiii;
          for(j=0;j<4;j++)
          {
            for(jj=j+1;jj<4;jj++)
            {
              if(is2[j]!=is2[jj])score+=gapbase[j]*gapbase[jj];
              score+=gapbase[j]*gapbase[jj];
            }
          }

          //printf("%d\n",score );
          if(score<min){min=score;is[0]=i;is[1]=ii;is[2]=iii;is[3]=iiii;}       
        }
      }
    }
  }

  //printf("%d\n\n",min );
  for(ii=0;ii<4;ii++)
  {
    for(i=0;i<Tiefe;i++)
    {
      if(is[ii])
      {
        if((*(Present_Column->Bases+i)==4 && *(Next_Column->Bases+i)==ii))
        {
          for(iii=0;iii<6;iii++) //subtracting W_Con
          { 
            if(iii!= *(Present_Column->Bases+i))Present_Column->w_con[iii]--; 
            if(iii!= *(Next_Column->Bases+i))Next_Column->w_con[iii]--; 
          }
          *(Next_Column->Bases+i)=4;
          *(Present_Column->Bases+i)=ii;
          for(iii=0;iii<6;iii++)  //adding W_Con
          { 
            if(iii!= *(Present_Column->Bases+i))Present_Column->w_con[iii]++; 
            if(iii!= *(Next_Column->Bases+i))Next_Column->w_con[iii]++; 
          }          
        }
      }
      else
      {
        if((*(Present_Column->Bases+i)==ii && *(Next_Column->Bases+i)==4))
        {
          for(iii=0;iii<6;iii++) //subtracting W_Con
          { 
            if(iii!= *(Present_Column->Bases+i))Present_Column->w_con[iii]--; 
            if(iii!= *(Next_Column->Bases+i))Next_Column->w_con[iii]--; 
          }          
          *(Next_Column->Bases+i)=ii;
          *(Present_Column->Bases+i)=4;
          for(iii=0;iii<6;iii++)  //adding W_Con
          { 
            if(iii!= *(Present_Column->Bases+i))Present_Column->w_con[iii]++; 
            if(iii!= *(Next_Column->Bases+i))Next_Column->w_con[iii]++; 
          }            
        }
      }
    }
  }
}

void Columns_Downdater(int reihe)
{
  #ifdef CURRENTFUNCTION
  printf("In Columns_Downdater\n");
  fflush(stdout);
  #endif 

  int j;
  char oldbase;
  Column* Present_Column=Null_Column;
  for(j=0;j<Breite;j++)
  {
    oldbase=*(Present_Column->Bases+reihe);

    int i;
    if(oldbase!=5)
    {
      for(i=0;i<6;i++) 
      {
        if(i!=oldbase)(*(Present_Column->w_con+i))-=1;      //Nur das hier kann zu unsigned <-> int overflow führen. Passiert aber nicht.
      }
    }
    Present_Column=Present_Column->Next;
  }

  #ifdef CURRENTFUNCTION
  printf("Out Columns_Downdater\n");
  fflush(stdout);
  #endif 
}

void Columns_Base_Downdater(int reihe)
{
  #ifdef CURRENTFUNCTION
  printf("In Columns_Base_Downdater\n");
  fflush(stdout);
  #endif 
  int j;
  Column* Present_Column=Null_Column;
  for(j=0;j<Breite;j++)
  {
    *(Present_Column->Bases+reihe)=5;  
    Present_Column=Present_Column->Next;
  }
  #ifdef CURRENTFUNCTION
  printf("Out Columns_Base_Downdater\n");
  fflush(stdout);
  #endif 
}

void Column_Updater(Column* Present_Column,char newbase,int reihe)
{ 
  #ifdef CURRENTFUNCTION
  printf("In Columns_Updater\n");
  fflush(stdout);
  #endif  

  int i;
  if(newbase!=5)
  {
    for(i=0;i<6;i++)
    {
      if(i!=newbase)(*(Present_Column->w_con+i))+=1;  
    }
  }
  *(Present_Column->Bases+reihe)=newbase;

  #ifdef CURRENTFUNCTION
  printf("Out Columns_Updater\n");
  fflush(stdout);
  #endif
}

void Column_Adder(Column* PreviousColumn,char newbase,int reihe)
{ 
  #ifdef CURRENTFUNCTION
  printf("In Column_Adder\n");
  fflush(stdout);
  #endif

  Column* New_Column;

  if(reserve<1)
  {
    New_Column=Guarded_Malloc(sizeof(Column));
    New_Column->Bases=Guarded_Malloc(sizeof(char)*Max_Seq_Anzahl);
    New_Column->w_con=Guarded_Malloc(sizeof(unsigned long)*6);
  }

  else
  {
    reserve--;
    New_Column=Reservoir[reserve];
  }


  Breite++;

  New_Column->Next=PreviousColumn->Next;
  New_Column->Previous=PreviousColumn;

  if (PreviousColumn!=Last_Column){(New_Column->Next)->Previous=New_Column;}
  (New_Column->Previous)->Next=New_Column;

/*  New_Column->Bases=Guarded_Malloc(sizeof(char)*Max_Seq_Anzahl);
  New_Column->w_con=Guarded_Malloc(sizeof(unsigned int)*6);*/

  unsigned int CovGapCount=0;
  unsigned int AlGapCount=0;
  int i;


  for(i=0;i<6;i++)*(New_Column->w_con+i)=0;

  if(PreviousColumn==Last_Column)
  {
    for(i=0;i<Tiefe;i++)  //entering the bases
    {
      if(i!=reihe)
      {        
        *(New_Column->Bases+i)=5;
        CovGapCount++;
      }  
    }
    Last_Column=New_Column;
  }

  else
  {
    for(i=0;i<Tiefe;i++)  //Entering the bases
    {
      if(i!=reihe)
      {  
        if( *((New_Column->Next)->Bases+i)==5 || *((New_Column->Previous)->Bases+i)==5 )
        {
          *(New_Column->Bases+i)=5;
          CovGapCount++;
        }
        else
        {
          *(New_Column->Bases+i)=4;
          AlGapCount++;
        }
      }  
    }
  }
  *(New_Column->Bases+reihe)=newbase;

  for(i=0;i<6;i++)
  {
    if(i!=newbase)(*(New_Column->w_con+i))+=1;
    if(i!=4)(*(New_Column->w_con+i))+=AlGapCount;
    //if(i!=5)(*(New_Column->w_con+i))+=CovGapCount;
  }
  //printf("NewColumn 5: %lu\n",*(New_Column->w_con+5) );

  #ifdef CURRENTFUNCTION
  printf("Out Column_Adder\n");
  fflush(stdout);
  #endif
}

int Backtracker(int length, int k)
{
  #ifdef CURRENTFUNCTION
  printf("In Backtracker\n");
  fflush(stdout);
  #endif

  int x,y,anf,u;
  int gapcount=1;

  Column* Present_Column;
  x=length-1;
  y=Breite-1;
  Present_Column=Last_Column;
  //printf("x=%d and y=%d\n",x,y);
  fflush(stdout);
  unsigned long oldscorey5;  

  int wayin=y;
  unsigned long best=MatrixOut(length-1,Breite-1);
  while( y>intmax(-1,Way[x]-bandwidthhalf) )  //y>-1  and sequences would be pushed out of the msa
  {
    if(MatrixOut(x,y)<best){best=MatrixOut(x,y);wayin=y;}
    y--;
  }

  y=wayin;
  //printf("wayin %d - %lu - %lu - %lu\n",wayin,MatrixOut(x,wayin+1),MatrixOut(x,wayin) ,MatrixOut(x,wayin-1));
  wayin=Breite-1;
  while(wayin>y)
  {
    Column_Updater(Present_Column,5,k);
    Present_Column=Present_Column->Previous;  
    wayin--;
  }


  while(x>-1 && y>-1)
  {
    anf=intmax(0,Way[x]-bandwidthhalf);

    if(MatrixOut(x,y)==MatrixOut(x,y-1)+Score(y,4))  //y-- has to come first
    {
      oldscorey5=Score(y,5); //The old score
      if(x==length-1){Column_Updater(Present_Column,5,k);}
      else
      {Column_Updater(Present_Column,4,k);}
      //if(y<Breite-1 && *(Present_Column->Next->Bases+k)==5){printf("BacktrackGapCrap5\n");exit(1);}
      Present_Column=Present_Column->Previous;      
      y--;
    }
    //Das hier neu jenseits der Seq kein Score mehr.
    else if(MatrixOut(x,y)==MatrixOut(x,y-1) && x==length-1)  //y-- has to come first
    {
      oldscorey5=Score(y,5); //The old score for diagnostics
      Column_Updater(Present_Column,5,k);
      Present_Column=Present_Column->Previous;      
      y--;
    }

    else if(MatrixOut(x,y)==MatrixOut(x-1,y-1)+Score(y,Seq_Bases[x]))
    {
      oldscorey5=Score(y,5); //alt Score(y+1,5) before update
      Column_Updater(Present_Column,Seq_Bases[x],k);
      Present_Column=Present_Column->Previous;
      if(Gaps[gapcount]==x && gapcount>0)gapcount--;       
      x--;
      y--;
    }

    else if(y>0 && MatrixOut(x,y)==MatrixOut(x-1,y)+u_intmax(Score(y,5),Score(y-1,5)))
    {
      //printf("y %d/%d Scorey %lu Score y-1 %lu\n",y,Breite,Score(y,5),Score(y-1,5) );
      Column_Adder(Present_Column,Seq_Bases[x],k);
      if(Gaps[gapcount]==x && gapcount>0)gapcount--;
      x--;
    }

    else
    {
      printf("\nStuff gone wrong\n");
      printf("seq %d, Breite %d\n",k,Breite );
      printf("x %d y %d length %d, oldscorey5 %lu\n",x,y,length,oldscorey5);
      printf("Mx-1y+Sc %lu %lu = %lu \n",MatrixOut(x-1,y),Score(y,5),MatrixOut(x-1,y)+Score(y,5) );
      printf("Mx-1y-1+Sc %lu %lu = %lu \n",MatrixOut(x-1,y-1),Score(y,Seq_Bases[x]),MatrixOut(x-1,y-1)+Score(y,Seq_Bases[x]) );
      printf("Mxy-1+Sc %lu %lu = %lu \n",MatrixOut(x,y-1),Score(y,4),MatrixOut(x,y-1)+Score(y,4));
      for(u=0;u<6;u++)printf("%lu ",w_con_array[y][u] );
      printf("\n");
      for(u=0;u<6;u++)printf("%lu ",w_con_array[y+1][u]);
      printf("\n");
      printf("Mxy%lu\n",MatrixOut(x,y));
      printf("wayin = %d, best = %lu\n",wayin,best );
      exit(1);
    }

/*    printf("\n");
    printf("Mx-1y+Sc %lu %lu \n",MatrixOut(x-1,y),Score(y,5) );
    printf("Mx-1y-1+Sc %lu %lu\n",MatrixOut(x-1,y-1),Score(y,Seq_Bases[x]) );
    printf("Mxy-1+Sc %lu %lu\n",MatrixOut(x,y-1),Score(y,4));    
    printf("M(%d)(%d)%lu, W %d\n",x,y,MatrixOut(x,y),Way[x]);*/
    if(MatrixOut(x,y)>9046744073709551615)exit(0);
  }
  
  while(y>-1)
  {
    Column_Updater(Present_Column,5,k);
    //if(y<Breite-1 && *(Present_Column->Next->Bases+k)==4){printf("BacktrackGapCrap4-%d-%d\n",x,y);exit(1);}
    Present_Column=Present_Column->Previous;     
    y--;
  }
  #ifdef PRINT
  printf("x=%d and y=%d\n",x,y);
  #endif

  #ifdef CURRENTFUNCTION
  printf("Out Backtracker\n");
  fflush(stdout);
  #endif

  return best;
}

void MatrixPrint(int length)
{
  int x,y;
  for(x=0;x<length;x++)
  {
    for(y=0;y<bandwidth;y++)
    {
      printf("%lu ",Matrix[x][y] );
    }
    printf("\n\n");
  }
}

unsigned int Matrix_Filler(int k)
{
  #ifdef CURRENTFUNCTION
  printf("In Matrix_Filler\n");
  fflush(stdout);
  #endif

  int x,y,anf,end;
  unsigned long eintrag;
  W_Con();  //W_Con first allows empty columns
  length=TheWay(k);  //The Way comes before Downdater and W_Con, otherwise W_Con bases are deleted
  Columns_Downdater(k);
  Columns_Base_Downdater(k);  //Only after TheWay bases are deleted
  

  #ifdef PRINT
  printf("Länge der Sequenz: %d\n",length);
  #endif

  if(length==0)return 0;

  //int gapcount=0;

  //Matrixinitialising is intrinsic in MatrixOut
  for(x=0;x<length;x++)
  {

  	anf=intmax(0,Way[x]-bandwidthhalf);
    end=intmin(Breite,anf+bandwidth);

    for(y=anf;y<end;y++)
    {
      //printf("%d %d\n",x,y);fflush(stdout);

      eintrag=MatrixOut(x-1,y-1)+Score(y,Seq_Bases[x]);
      eintrag=u_intmin(eintrag,MatrixOut(x,y-1)+Score(y,4));
      if(y>0 && y<Breite-1)     //No Column_Add at the end.
      {
        eintrag=u_intmin(eintrag,MatrixOut(x-1,y)+u_intmax(Score(y,5),Score(y-1,5)) );
      }  //Score(y,5) number of not ' '. Only an approximation of y-0.5.    
      
      MatrixIn(x,y,eintrag);  

    }
  }
  #ifdef PRINT
  printf("MatrixFilling done.\n");
  fflush(stdout); 
  printf("This is the Score: %u\n",MatrixOut(length-1,Breite-1));
  fflush(stdout); 
  #endif

  //MatrixPrint(length);
  #ifdef CURRENTFUNCTION
  printf("Out Matrix_Filler\n");
  fflush(stdout);
  #endif

  return Backtracker(length,k);

  //return MatrixOut(length-1,Breite-1);

}



void Freer()
{
  #ifdef PRINT
  printf("Freeing Stuff.\n");
  fflush(stdout);
  #endif

  int i;
  Column* Present_Column;
  Present_Column=Null_Column;
  for(i=0;i<Breite-1;i++)
  {
    free(Present_Column->Bases);
    free(Present_Column->w_con);
    Present_Column=Present_Column->Next;
    free(Present_Column->Previous);
  }
  free(Present_Column);
  free(w_con_array);
}

void MMA_Auslesen(char *outputfile)
{
  Chars[0]='A';
  Chars[1]='C';
  Chars[2]='G';
  Chars[3]='T';
  Chars[4]='-';
  Chars[5]=' ';
  int i,j;
  FILE *datei;
  datei=fopen(outputfile,"w");

  if(NULL == datei)
  {
    printf("DateiVerbratei!\n");
    exit(1);
  }

  #ifdef PRINT
  printf("File opened\n");
  printf("%dx%d MMA\n",Tiefe,Breite);
  fflush(stdout);
  #endif

  Column* Present_Column;
  for(j=0;j<Tiefe;j++)
  {
    Present_Column=Null_Column;
    for(i=0;i<Breite;i++)
    {
      fprintf(datei,"%c", Chars[(int)*(Present_Column->Bases+j)] );
      //printf("%c", Chars[*(Present_Column->Bases+j)] );
      Present_Column=Present_Column->Next;
    }
    fprintf(datei,"\n"); 
  }
  fclose(datei); 

  #ifdef PRINT
  printf("File filled\n");
  fflush(stdout);  
  #endif
}

void Help()
{
  printf("Usage: ./PW_ReAligner MApath\n");
  printf("Flags:\n");
  printf("-o msa_path    Path of the refined multiple sequence alignment. Default: SimulatedMSAreal.\n");
  printf("-b <1000>      The width of the band that is calculated in the alignment matrix.\n");
  exit(0);
}


int main(int argc, char *argv[])
{


  char *MApath_p;
  if(argc<2){printf("Usage: ./PW_ReAligner MApath\n");exit(0);}
  MApath_p=argv[1];

  //Specifying the output path
  char outputfile[]="MSAreal";
  char *output_p=&outputfile[0];
 

  int i;

  bandwidth=1000;
  bandwidthhalf=bandwidth/2;

  for(i=1;i<argc;i++)
  {

    if(argv[i][0]=='-' && argv[i][1]=='o')
    {
      printf("%s\n",argv[i]);
      if(i+1<argc)output_p=argv[i+1];
    }  

    if(argv[i][0]=='-' && argv[i][1]=='b')
    {
      bandwidth=atoi(argv[i+1]);
      bandwidthhalf=bandwidth/2;
    }

    if(argv[i][0]=='-' && argv[i][1]=='h')
    {
      Help();
    } 
  }

  printf("output file: %s\n",output_p );
  printf("bandwidth %d\n",bandwidth );


  MMA_Einlesen(MApath_p);

  EntAlGapper();
  
  printf("Rows %d, Columns %d.\n",Tiefe,Breite);
     
  W_Con();
  W_Con_Checker();


  
  BestMille=-1;
  OverallScorePrint();


  #ifdef PRINT
  printf("Breadth: %d\n",Breite);
  #endif
  //printf("ReihenScore: %d\n",ReihenScore(11));


  int k;
  unsigned long newscore,oldscore,reihenscore;
  i=0; 
  int nichtverbessert=0;

  clock_t begin = clock();
  
  while(i<10000)
  {
    if(0) //i==5) //bandwidth reduction
    {
      bandwidth-=((bandwidth*4)/5);
      bandwidthhalf=bandwidth/2;
    }
    if(0)  //i==15) //bandwidth reduction
    {
      bandwidth/=2;
      bandwidthhalf=bandwidth/2;
    }    

    
    for(k=0;k<Tiefe;k++) 
    {
      int oldbreite=Breite;
      //EntAlGapper();

      //oldscore=ReihenScore(k);

      newscore=Matrix_Filler(k);

      //printf("%d\n",(int)((newscore*100)/(Lengths[k]*Tiefe)) );fflush(stdout);
      //reihenscore=ReihenScore(k);
      //printf("%d:%d %lu->%lu(%lu),l%d, Breite %d\n",i,k,oldscore,newscore,reihenscore,TheWay(k),Breite);

      if(0) //(newscore>oldscore || reihenscore>oldscore)
      {
        printf("Error in round %d Alignment %d\n",i,k);
        printf("Breadth %d -> %d\n",oldbreite,Breite );
        printf("%lu %lu %lu\n",oldscore,newscore,reihenscore);
        printf("%d ",Matrix_Filler(k));
        printf("%lu ",ReihenScore(k));
        printf("%d \n",Matrix_Filler(k));
        printf("length %d\n",TheWay(k) );
        printf("Way %d -> %d\n",Way[0],TheWay(k)-1 );
        printf("%d -> %d/%d\n",Way[0],Way[Lengths[k]-1],Breite);
        exit(0);
      }

      
      if(k%500==0)
      { 
        W_Con();
        #ifdef PRINT
        printf("W_Con_Checker:%d\n",W_Con_Checker());
        printf("Breite: %d\n",Breite);
        fflush(stdout);
        #endif
      }
      #ifdef PRINT
      printf("Round %d, row %d \n",i,k);
      #endif
      //printf("Durchgang %d, Reihe %d \n",i,k);
   
    }

    //EntAlGapper();

    if(OverallScorePrint()){MMA_Auslesen(output_p);nichtverbessert=0;}
    else{nichtverbessert+=1;if(nichtverbessert>0)i=100000;}
    //printf("%lu\n",OverallScore() );

    i++;

  }

  clock_t end_t = clock();
  double total_t = ((double)(end_t - begin) / CLOCKS_PER_SEC)/60;
  printf("Total time: %f min.\n", total_t  );

  EntAlGapper();
  if(OverallScorePrint())MMA_Auslesen(output_p);

  Freer();

  exit(0);
}













