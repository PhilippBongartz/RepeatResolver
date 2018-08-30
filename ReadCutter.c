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

#define maxlength2 35000   // max read length
#define maxlength1 40000   // max template piece length


//#define PRINT

/*
This tool divides the reads into repeat sections and unique sections. 
*/




void *Guarded_Malloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL)
    { fprintf(stderr,"\nError in LargeScaleVars\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}

void *Guarded_Realloc(void *p, size_t size)
{ p = realloc(p,size);
  if (p == NULL)
    { fprintf(stderr,"\nError in LargeScaleVars:\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}

int intmin(int x, int y){if(x<y){return x;} return y;}

int intmax(int x, int y){if(x>y){return x;} return y;}

int intabs(int a)
{
  if(a<0){return -a;}
  return a;
}

//ReadingFasta -> Reads: With offset to make it possible to read in reads by readnumber:
off_t offset=0; //global
int readcount=-1;
char *Read;
int readlength=0;
char *buffer;
int readanzahl;
int ReadingFasta(char *Path)
{
  char *s;
  FILE * File=NULL;
  size_t len=70000;

  File=fopen(Path,"r");
  if(fseek(File,offset,SEEK_SET)){printf("Problem.\n");exit(1);}

  int i;

  //printf("%p\n",File );
  fflush(stdout);

  if(File==NULL){exit(1);}

  //File+=offset;
  int bytes=0;
  int basecount=0;
  int nextread=0;
  while(nextread<2)
  {
    s = fgets(buffer, len, File);
    if(s==NULL){offset=0;return 1;}

    if(buffer[0]=='>')
    {
      nextread+=1;
      i=0;
      while(buffer[i]!='\n' && nextread<2)  //the start of the read
      {
        i+=1;
      }
      bytes+=i+1;  //+1 because of \n
    }
    else  //the bases are read in
    {
      //printf("\n\n%d\n",readcount);
      i=0;
      while(buffer[i]!='\n')
      {
        if(buffer[i]=='A' || buffer[i]=='a')Read[basecount]='a';
        else if(buffer[i]=='C' || buffer[i]=='c')Read[basecount]='c';
        else if(buffer[i]=='G' || buffer[i]=='g')Read[basecount]='g';
        else if(buffer[i]=='T' || buffer[i]=='t')Read[basecount]='t';
        else {basecount--;}        
        basecount++;
        //printf("%c ",buffer[i]);
        i++;
      }
      bytes+=i+1;  //+1 because of \n
    }

  }
  readcount++;

  fclose(File);
  offset+=bytes-1;  //-1 to set it back to '>'. 
  //printf("offset: %lld\n",offset );
  Read[basecount]='\n';

  readlength=basecount;

  #ifdef PRINT
  printf("Readlength %d\n",basecount );
  printf("Read %d\n",readcount );
  fflush(stdout);
  #endif
  
  return 0;
}


char BaseInverse[300];
void ReadInversion()
{
  int i;
  char temp;
  for(i=readlength/2;i>-1;i--)
  {
    temp=Read[i];
    Read[i]=BaseInverse[(int)Read[readlength-i-1]];
    Read[readlength-i-1]=BaseInverse[(int)temp];
  }
}

//This reads in the template
char *Template;
int templatelength=0;
void ReadingTemplate(char *Path)
{
  char *s;
  FILE * File=NULL;
  size_t len=70000;
  File=fopen(Path,"r");

  int i;

  if(File==NULL){printf("No template.\n");fflush(stdout);exit(1);}

  int basecount=0;
  //int nextread=0;

  while((s = fgets(buffer, len, File)) != NULL)
  {
    if(buffer[0]!='>')
    {
      i=0;
      while(buffer[i]!='\n')
      {
        if(buffer[i]=='A' || buffer[i]=='a')Template[basecount]='a';
        else if(buffer[i]=='C' || buffer[i]=='c')Template[basecount]='c';
        else if(buffer[i]=='G' || buffer[i]=='g')Template[basecount]='g';
        else if(buffer[i]=='T' || buffer[i]=='t')Template[basecount]='t';
        else {basecount--;}
        basecount++;
        i++;
      }
    }
  }

  fclose(File);
  //Template[basecount]='\n';
  templatelength=basecount;

  //wrap around
  for(i=0;i<templatelength;i++){Template[basecount+i]=Template[basecount];}
}

// Speed up possibility for a future version, not currently used. 
int Bitvectoraligner(char *shortstring, char *longstring)
{
  int length1=strlen(shortstring);
  int length2=strlen(longstring);

  long unsigned int B[4];
  long unsigned int mask;
  int i;
  int score[maxlength2];
  //Preprocessing:
  B[0]=0;
  B[1]=0;
  B[2]=0;
  B[3]=0;

  for(i=0;i<length1;i++)
  {
  	mask=1;
  	mask=mask<<i;
  	if(shortstring[i]=='a' || shortstring[i]=='A')
  	{
  	  B[0]=B[0] | mask;
  	}
  	if(shortstring[i]=='c' || shortstring[i]=='C')
  	{
  	  B[1]=B[1] | mask;  		
  	}
  	if(shortstring[i]=='g' || shortstring[i]=='G')
  	{
  	  B[2]=B[2] | mask;  		
  	}
  	if(shortstring[i]=='t' || shortstring[i]=='T')
  	{
  	  B[3]=B[3] | mask;  		
  	}
  }
  score[0]=length1;

  //Alignment:
  long unsigned int X;
  long unsigned int D[maxlength2];
  long unsigned int HP[maxlength2];
  long unsigned int HN[maxlength2];
  long unsigned int VP[maxlength2];
  long unsigned int VN[maxlength2];

  VP[0]=0;
  VP[0]=~VP[0];

  VN[0]=0;

  for(i=1;i<length2+1;i++)
  {
  	if(longstring[i-1]=='a' || longstring[i-1]=='A')X=B[0] | VN[i-1];
  	if(longstring[i-1]=='c' || longstring[i-1]=='C')X=B[1] | VN[i-1];
  	if(longstring[i-1]=='g' || longstring[i-1]=='G')X=B[2] | VN[i-1];
  	if(longstring[i-1]=='t' || longstring[i-1]=='T')X=B[3] | VN[i-1];

  	D[i]=((VP[i-1]+(X&VP[i-1]))^VP[i-1])|X;

  	HN[i]=VP[i-1]&D[i];
  	HP[i]=VN[i-1]| ~(VP[i-1]|D[i]);
    mask=1;
  	X=HP[i]<<1 | mask;
  	VN[i]=X&D[i];
  	VP[i]=(HN[i]<<1)|~(X|D[i]);

  	mask=1;
  	mask=mask<<(length1-1);

  	if(HP[i]&mask)
    {
      score[i]=score[i-1]+1;
    }
  	else
    {
      if(HN[i]&mask)score[i]=score[i-1]-1;
      else score[i]=score[i-1];
    }
  }

  //Backtracking:
  int eslength=0;
  char es[maxlength2*2];
  int j;
  j=length1-1;
  i=length2;
  while(i>0 && j>-1)
  {
    //printf("i:%d j:%d\n",i,j);
    mask=1;
    mask=mask<<j;
    if(D[i]&mask && shortstring[j]==longstring[i-1])
    {
      i--;
      j--;
      es[eslength]='m';
      eslength++;
    }
    else if(HP[i]&mask)
    {
      i--;
      es[eslength]='d';
      eslength++;
    }
    else if(VP[i]&mask)
    {
      j--;
      es[eslength]='i';
      eslength++;
    }
    else
    {
      j--;
      i--;
      es[eslength]='s';
      eslength++;
    }    
  }
  printf("j %d, i %d\n",j,i);
  while(i>0)
  {
     i--;
     es[eslength]='d';
     eslength++;  	
  }
  printf("Editscript:");
  for(i=0;i<eslength;i++)printf("%c",es[i]);
  printf("\n");  

	//Adapting the output to ./Aligner:
  int alscore=0;

  for(i=0;i<eslength;i++)
  {
    if(es[i]!='m')alscore++;
  }
  alscore-=(length2-length1);

  printf("%d\n",alscore);
  printf("\n");

  for(i=eslength-1;i>=0;i--)
  {
    if(es[i]=='m')printf("m");
    if(es[i]=='s')printf("s");
    if(es[i]=='i')printf("d");
    if(es[i]=='d')printf("i");
  }
  printf("\n");


  return alscore;
}

// Implementation of a variation on Needleman-Wunsch
int **Matrix;    //The matrix without the left column and upper row, those can be accessed with -1 for intitialisation of edge cases.
int **Schatten1;   //The whole alignment Matrix
int **Schatten2;   //The matrix except the left column
int IntoAligner(char *shortstring, char *longstring, int length1, int length2)
{
  //int length1=strlen(shortstring);
  //int length2=strlen(longstring);

  int x,y;


  //Initialisation:
  for(x=-1;x<length1;x++)Matrix[x][-1]=x+1;
  for(y=0;y<length2;y++)Matrix[-1][y]=0;    //Seq1 is aligned into Seq2.

  //Filling the matrix 
  int m;
  for(x=0;x<length1;x++)
  {
  	for(y=0;y<length2;y++)
  	{
  	  m=1;	
  	  if(shortstring[x]==longstring[y])m=0;
  	  Matrix[x][y]=Matrix[x-1][y-1]+m;
      Matrix[x][y]=intmin(Matrix[x-1][y]+1,Matrix[x][y]);

      if(0) //x==length1-1)  
        {Matrix[x][y]=intmin(Matrix[x][y-1],Matrix[x][y]);}
      else 
        {Matrix[x][y]=intmin(Matrix[x][y-1]+1,Matrix[x][y]);}
  	  
  	}
  }

  #ifdef PRINT
  printf("%d\n",Matrix[length1-1][length2-1]);
  fflush(stdout);
  #endif

  //Backtracking
  char *EditScript;
  int count=0;
  EditScript=Guarded_Malloc(sizeof(char)*(length1+length2));
  x=length1-1;
  y=length2-1;

  int min,i,einstieg_y;
  einstieg_y=y;
  min=Matrix[x][einstieg_y];
  for(i=length2-1;i>0;i--)
  {
  	if(Matrix[x][i]<min)
  	{
  		min=Matrix[x][i];
  		einstieg_y=i;
  	}
  }

  #ifdef PRINT
  printf("Einstieg %d\n",einstieg_y );
  fflush(stdout);
  #endif

  y=einstieg_y;

  while(x>-1 && y>-1)
  {
  	m=1;	
  	if(shortstring[x]==longstring[y])m=0;

  	if(Matrix[x][y]==Matrix[x-1][y-1]+m)
  	{
  	  if(m) { EditScript[count]='s';}
  	  else  {EditScript[count]='m';}
  	  count++;
  	  x--;
  	  y--;
  	}
    else if(Matrix[x][y]==Matrix[x][y-1] && x==length1-1)
    {
      EditScript[count]='i';  
      count++;  
      y--;
    }    
  	else if(Matrix[x][y]==Matrix[x][y-1]+1 && x!=length1-1)
  	{
  	  EditScript[count]='i';	
  	  count++;	
  	  y--;
  	}
  	else if(
  		Matrix[x][y]==Matrix[x-1][y]+1)
  	{
  	  EditScript[count]='d';	
  	  count++;
  	  x--;  		
  	}
  	else
  	{
  	  printf("Backtrackingproblem.\n");
  	  exit(1);
  	}	
  	//printf("x:%d y:%d Matr:%d\n",x,y,Matrix[x][y]);
  }

  while(x>-1)
  {
  	EditScript[count]='d';	
    count++;
	x--;  	  	
  } 
  while(y>-1) 	
  {
  	EditScript[count]='i';	
  	count++;	
  	y--;  	
  }

  //inverting the edit script:
  for(x=0;x<count/2;x++)
  {
  	m=EditScript[x];
  	EditScript[x]=EditScript[count-x-1];
  	EditScript[count-x-1]=m;
  }

  #ifdef PRINT
  for(x=0;x<count;x++)printf("%c",EditScript[x]);
  printf("\n"); 
  fflush(stdout);
  #endif


  return min;
}

//This function detects all occurences of a short string in a longer string:
int Positions[100];
int pos_count;
void Occurrence(char *shortstring, char *longstring, int length1, int length2, int score_cutoff)
{
  //int length1=strlen(shortstring);
  //int length2=strlen(longstring);

  int x,y;


  //Initialisiation:
  for(x=-1;x<length1;x++)Matrix[x][-1]=x+1;
  for(y=0;y<length2;y++)Matrix[-1][y]=0;    //Seq1 is aligned into Seq2.

  //Filling the matrix
  int m;
  for(x=0;x<length1;x++)
  {
    for(y=0;y<length2;y++)
    {
      m=1;  
      if(shortstring[x]==longstring[y])m=0;
      Matrix[x][y]=Matrix[x-1][y-1]+m;
      Matrix[x][y]=intmin(Matrix[x-1][y]+1,Matrix[x][y]);

      if(0) //x==length1-1)    
        {Matrix[x][y]=intmin(Matrix[x][y-1],Matrix[x][y]);}
      else 
        {Matrix[x][y]=intmin(Matrix[x][y-1]+1,Matrix[x][y]);}
      
    }
  }

  //printf("%d\n",Matrix[length1-1][length2-1]);

  // Each detected occurence is written into Positions.
  int on=0;
  int min,i,einstieg_y,lastmin;
  lastmin=100000;
  min=100000;
  pos_count=0;
  x=length1-1;

  for(i=length2-1;i>0;i--)
  {
    if(Matrix[x][i]<score_cutoff)
    {
      on=1;  //We have found an occurence, now we look for the minimal score of this occurence
    }
    else
    {
      if(on)
      {
        if(pos_count>0 && Positions[pos_count-1]-einstieg_y>length1/2) // distance to the last occurence is big enough: this is a new one
        {
          Positions[pos_count]=einstieg_y;
          pos_count++;
        }
        else if(pos_count>0 && Positions[pos_count-1]-einstieg_y<=length1/2) // still the same occurence
        {
          if(lastmin>min){Positions[pos_count-1]=einstieg_y;}  // better score -> replaced
        }
        else if(pos_count==0)
        {
          Positions[pos_count]=einstieg_y;
          pos_count++;
        }
      }
      on=0;
      lastmin=min;
      min=100000;
    }

    if(on && Matrix[x][i]<min)
    {
      min=Matrix[x][i];
      einstieg_y=i;
    }
  }
}



// Detecting all parts of the template in a read: 
int *Part_Indices;
int *Part_Positions;
int **OrderMatrix;
int ***DistanceHistos;
int **CuttingPoints;
int *Cutting_Number;
int *Inverted;
int allpartscount;
void FullAnalysis(int parts, int overlap, double error_cutoff, int wiggleroom)
{
  int steps=templatelength/parts;
  int len=steps+overlap;
  int cutoff=(int)(((double)len)*error_cutoff);


  #ifdef PRINT
  printf("steps %d, len %d, cutoff %d\n",steps,len,cutoff );
  fflush(stdout);
  #endif


  int i,j,k,temp;

  //Mapping of all parts into the reads
  k=0;
  for(i=0;i<parts;i++)
  {
    if(i==0 || i==parts-1)  // Otherwise runtime is prohibitive:
    {Occurrence(&Template[i*steps],&Read[0],len,readlength, cutoff);}
    //printf("Positions: ");
    for(j=0;j<pos_count;j++)
    {
      Part_Positions[k]=Positions[j];
      Part_Indices[k]=i;
      k++;
      //printf("%d ", Positions[j]);
    }
    //printf("\n");    
  }
  allpartscount=k;

  //Sorting the mappings by position
  for(i=0;i<k;i++)
  {
    for(j=0;j<k;j++)
    {
      if(Part_Positions[i]<Part_Positions[j])
      {
        temp=Part_Positions[i];
        Part_Positions[i]=Part_Positions[j];
        Part_Positions[j]=temp;

        temp=Part_Indices[i];
        Part_Indices[i]=Part_Indices[j];
        Part_Indices[j]=temp;
      }
    }
  }

  #ifdef PRINT
  for(i=0;i<k;i++)printf("Part %d at %d\n",Part_Indices[i],Part_Positions[i]);
  fflush(stdout);
  #endif


  //Filling the order matrix:
  for(i=0;i<k-1;i++)
  {
    if( abs(Part_Positions[i]+len-Part_Positions[i+1])<wiggleroom)
    {
      OrderMatrix[Part_Indices[i]][Part_Indices[i+1]]++;
    }
  }

  //Filling the DistanceHistos:

  for(i=0;i<k-1;i++)
  {
    j=(int)(log(1+abs(Part_Positions[i]+len-Part_Positions[i+1]))/log(1.5));
    //printf("%d %d -> %d = %d\n",Part_Positions[i],Part_Positions[i+1],abs(Part_Positions[i]+len-Part_Positions[i+1]),j );
    DistanceHistos[Part_Indices[i]][Part_Indices[i+1]][j]++;
  }



  //Cutting the reads up: CuttingPoints, Cutting_Number
  Cutting_Number[readcount]=0;
  for(i=0;i<k;i++)
  {
    if(Part_Indices[i]==parts-1 && Part_Positions[i]>len && readlength-Part_Positions[i]>len)
    {
      CuttingPoints[readcount][Cutting_Number[readcount]]=Part_Positions[i];
      Cutting_Number[readcount]++;
    }
  }

  #ifdef PRINT
  printf("Cuts: ");
  for(i=0;i<Cutting_Number[readcount];i++)
  {
    printf("%d ",CuttingPoints[readcount][i] );
  }
  printf("\n");
  fflush(stdout);
  #endif



  if(parts>1)
  {
    // This is a try to calculate more robust cutting points:
    j=0;
    for(i=0;i<k;i++)
    {
      if(Part_Indices[i]==parts-1 && Part_Positions[i]>len && readlength-Part_Positions[i]>len)
      {
        CuttingPoints[readcount][j]=Part_Positions[i];
        j++;
      }
    }
    for(i=0;i<k;i++)
    {
      if(Part_Indices[i]==0 && Part_Positions[i]-len>len && readlength-(Part_Positions[i]-len)>len)
      {
        CuttingPoints[readcount][j]=Part_Positions[i]-len;
        j++;
      }
    }
    for(i=0;i<k;i++)
    {
      if(Part_Indices[i]==parts-2 && Part_Positions[i]+len>len && readlength-(Part_Positions[i]+len)>len)
      {
        CuttingPoints[readcount][j]=Part_Positions[i]+len;
        j++;
      }
    }
    for(i=0;i<k;i++)
    {
      if(Part_Indices[i]==1 && Part_Positions[i]-2*len>len && readlength-(Part_Positions[i]-2*len)>len)
      {
        CuttingPoints[readcount][j]=Part_Positions[i]-2*len;
        j++;
      }
    }

    //Searching for the best cutting points:
    Cutting_Number[readcount]=0;
    //First the first:
    for(i=0;i<j;i++)
    {
      if(CuttingPoints[readcount][i]<templatelength+templatelength/2)
      {
        CuttingPoints[readcount][Cutting_Number[readcount]]=CuttingPoints[readcount][i];
        Cutting_Number[readcount]++;
        i=j;
      }
    }
    //Then the rest:
    for(k=0;k<60;k++)
    {
      for(i=0;i<j;i++)
      {
        if(CuttingPoints[readcount][Cutting_Number[readcount]-1]+templatelength/2<CuttingPoints[readcount][i] && CuttingPoints[readcount][i]<CuttingPoints[readcount][Cutting_Number[readcount]-1]+templatelength+templatelength/2 ) 
        {
          CuttingPoints[readcount][Cutting_Number[readcount]]=CuttingPoints[readcount][i];
          Cutting_Number[readcount]++;
          i=j;
        }
      } 
    }

    #ifdef PRINT
    printf("Cuts: ");
    for(i=0;i<Cutting_Number[readcount];i++)
    {
      printf("%d ",CuttingPoints[readcount][i] );
    }
    printf("\n");
    fflush(stdout);
    #endif


  }

}


// Allocation of most of the memory
int loglen;
void MatrixMemory(int parts, int readanzahl)
{
  //Memory allocation
  Schatten1=Guarded_Malloc(sizeof(int*)*(maxlength1+1));
  int x;
  for(x=0;x<maxlength1+1;x++)*(Schatten1+x)=Guarded_Malloc(sizeof(int)*(maxlength2+1));
  Schatten2=Guarded_Malloc(sizeof(int*)*(maxlength1+1));
  for(x=0;x<maxlength1+1;x++)*(Schatten2+x)=(*(Schatten1+x)+1); 
  Matrix=Schatten2+1;  

  Part_Indices=Guarded_Malloc(sizeof(int)*1000);
  Part_Positions=Guarded_Malloc(sizeof(int)*1000);

  int i;
  OrderMatrix=Guarded_Malloc(sizeof(int*)*parts);
  for(x=0;x<parts;x++)
  {
    OrderMatrix[x]=Guarded_Malloc(sizeof(int)*parts);
    for(i=0;i<parts;i++)OrderMatrix[x][i]=0;
  }

  int j;
  loglen=1+(int)(log(100000)/log(1.5));

  #ifdef PRINT
  printf("loglen: %d\n",loglen );
  fflush(stdout);
  #endif


  DistanceHistos=Guarded_Malloc(sizeof(int**)*parts);
  for(i=0;i<parts;i++)
  {
    DistanceHistos[i]=Guarded_Malloc(sizeof(int*)*parts);
    for(j=0;j<parts;j++)
    {
      DistanceHistos[i][j]=Guarded_Malloc(sizeof(int)*loglen);
      for(x=0;x<loglen;x++)DistanceHistos[i][j][x]=0;
    }
  }

  CuttingPoints=Guarded_Malloc(sizeof(int*)*readanzahl);
  Cutting_Number=Guarded_Malloc(sizeof(int)*readanzahl);
  Inverted=Guarded_Malloc(sizeof(int)*readanzahl);

  for(i=0;i<readanzahl;i++)CuttingPoints[i]=Guarded_Malloc(sizeof(int)*60);
  for(i=0;i<readanzahl;i++)Cutting_Number[i]=0;
  for(i=0;i<readanzahl;i++)Inverted[i]=0;

  BaseInverse['a']='t';
  BaseInverse['t']='a';
  BaseInverse['c']='g';
  BaseInverse['g']='c';
  BaseInverse['A']='T';
  BaseInverse['T']='A';
  BaseInverse['C']='G';
  BaseInverse['G']='C';

  Template=Guarded_Malloc(sizeof(char)*70000);

}

void Print_OrderMatrix(int parts)
{
  int i,x;
  for(x=0;x<parts;x++)
  {
    for(i=0;i<parts;i++)
    {
      printf("%6d ",OrderMatrix[x][i]);
    }
    printf("\n");
  }

}

void Print_DistanceHistos(int parts, int relevant)
{
  int i,j,k,count;
  for(j=0;j<parts;j++)
  {
    for(i=0;i<parts;i++)
    {
      count=0;
      for(k=0;k<loglen;k++){count+=DistanceHistos[j][i][k];}

      if(count>relevant)
      {
        printf("%d -> %d: ", j,i);
        for(k=0;k<loglen;k++){printf("%d ",DistanceHistos[j][i][k]);}
        printf("\n");
      }
    }
  }  
}

int ReadCounter(char *Path)
{
  char *s;
  FILE * File=NULL;
  size_t len=70000;
  File=fopen(Path,"r");
  int count=0;
  if(File==NULL){printf("No Reads.\n");fflush(stdout);exit(1);}
  while((s = fgets(buffer, len, File)) != NULL)
  {
    if(buffer[0]=='>')count++;
  }
  fclose(File);
  return count;
}

void Help()
{
  printf("Usage: ./ReadCutter template.fasta reads.fasta\n");
  printf("Flags:\n");
  printf("-p <20>    determines the number of parts into which the template is cut and which are mapped into the reads.\n");
  printf("-e <0.30>  is the mapping error cutoff being used to detect occurences of parts in reads\n");
  printf("-w <150>   restricts how far apart and close together mappings can be.\n");
  printf("-l <0>     is used to create parts which overlap by l bases.\n");
  exit(0);
}

//This cuts reads into repeat instances: 
int seqcount;
void OutputOfCuts(FILE *Seqfasta)
{
  #ifdef PRINT
  printf("OutputOfCuts Read %d\n",readcount);
  fflush(stdout);
  #endif

  int i,j;
  j=0;
  fprintf(Seqfasta,">\n");

  for(i=0;i<readlength;i++)
  {
    if(i==CuttingPoints[readcount][j] && j<Cutting_Number[readcount])
    {
      #ifdef PRINT
      printf("%d \n",i );
      fflush(stdout);
      #endif

      fprintf(Seqfasta,"\n>\n");    
      j++;
    }
    fprintf(Seqfasta,"%c",Read[i]);
  }
  fprintf(Seqfasta,"\n");
}


//This function outputs the information which sequences are in which read
int seqcount;
void OutputOfReadSeqInfo(FILE *ReadSeqFile)
{

  #ifdef PRINT
  printf("OutputOfReadSeqInfo\n");
  fflush(stdout);
  #endif

  int i,j;
  seqcount=0;
  for(i=0;i<readanzahl;i++)
  {
    for(j=0;j<Cutting_Number[i]+1;j++)
    {
      fprintf(ReadSeqFile, "%d ",seqcount );
      seqcount++;
    }
    fprintf(ReadSeqFile, "\n");
  }
}

int main(int argc, char *argv[])
{

  char *Templatepath_p;
  char *Readspath_p;
  if(argc<2){printf("Usage: ./ReadCutter template.fasta Reads.fasta\n");exit(0);}
  Readspath_p=argv[2];
  Templatepath_p=argv[1];

  //Building a filename with the right prefix Type+'_05perc_30kb_Template.fasta
  char outputfile[300];
  char readseqfile[300];

  int i=0;
  char TemplateSuffix[]="Template.fasta";
  char outputsuffix[]="Seq.fasta";
  char readseqsuffix[]="ReadSeqInfo";
  char DataPrefix[300];

  while(strcmp(Templatepath_p+i, TemplateSuffix)!=0 && i<strlen(Templatepath_p))
  {
    DataPrefix[i]=Templatepath_p[i];
    i++;
  }
  if(strcmp(Templatepath_p+i, TemplateSuffix)!=0){DataPrefix[0]='\0';}
  DataPrefix[i]='\0';

  strcpy(outputfile,DataPrefix);
  strcat(outputfile,outputsuffix);

  strcpy(readseqfile,DataPrefix);
  strcat(readseqfile,readseqsuffix);

  printf("outputfile: %s\n",outputfile);
  printf("readseqfile: %s\n",readseqfile);
  //////////////////////////

  //the output paths: 
  //char outputfile[]="Seq.fasta";
  char *output_p=&outputfile[0];
 
  //char readseqfile[]="ReadSeqInfo";
  char *readseq_p=&readseqfile[0];

  //arguments: -o outputfile, -p parts into which the template is cut, -h help, , -e error_cutoff,
  // -w wiggleroom, -l overlap
  int parts=60;
  int overlap=0;
  double error_cutoff=0.30;
  int wiggleroom=150;

  // parameters and paths
  for(i=1;i<argc;i++)
  {

    if(argv[i][0]=='-' && argv[i][1]=='o')
    {
      //printf("%s\n",argv[i]);
      if(i+1<argc)output_p=argv[i+1];
    }  

    if(argv[i][0]=='-' && argv[i][1]=='r')
    {
      //printf("%s\n",argv[i]);
      if(i+1<argc)readseq_p=argv[i+1];
    }  

    if(argv[i][0]=='-' && argv[i][1]=='p')
    {
      parts=atoi(argv[i+1]);
    }    

    if(argv[i][0]=='-' && argv[i][1]=='l')
    {
      overlap=atoi(argv[i+1]);
    }   

    if(argv[i][0]=='-' && argv[i][1]=='w')
    {
      wiggleroom=atoi(argv[i+1]);
    }   

    if(argv[i][0]=='-' && argv[i][1]=='e')
    {
      error_cutoff=atof(argv[i+1]);
    }   

    if(argv[i][0]=='-' && argv[i][1]=='h')
    {
      Help();
    }    
  }


  //printf("outputfile: %s\n",output_p );
  printf("parts %d, overlap %d, wiggleroom %d, error_cutoff %f\n",parts,overlap,wiggleroom,error_cutoff);

  buffer=Guarded_Malloc(sizeof(char)*70002);
  Read=Guarded_Malloc(sizeof(char)*70000);

  readanzahl=ReadCounter(Readspath_p);

  printf("read count %d\n",readanzahl );fflush(stdout);

  //Memory allocation:
  MatrixMemory(parts,readanzahl);


  ReadingTemplate(Templatepath_p);
  printf("template length %d\n",templatelength );


  //ReadingFasta(Readspath_p);
  int prozent=5;
  for(i=0;i<readanzahl;i++)
  {
    ReadingFasta(Readspath_p);
    FullAnalysis(parts, overlap, error_cutoff, wiggleroom);
    if(0) //allpartscount<1)
    {
      #ifdef PRINT
      printf("ReadInversion\n");
      #endif
      Inverted[readcount]=1;
      ReadInversion();
      FullAnalysis(parts, overlap, error_cutoff, wiggleroom);
    }
    if(i%10==0)
    {
      #ifdef PRINT
      Print_OrderMatrix(parts);
      Print_DistanceHistos(parts, 10);
      #endif
    }

    if(i*100/readanzahl>prozent)
    {
      printf("%d %% done.\n",prozent );
      prozent+=5;
    }
  }

  //Output: 
  printf("Outputting results.\n");

  seqcount=0;
  FILE *Seqfasta;
  Seqfasta=fopen(output_p,"w");

  readcount=-1;
  offset=0;

  for(i=0;i<readanzahl;i++)
  {
    ReadingFasta(Readspath_p);
    if(Inverted[readcount])
    {
      ReadInversion();
    }
    OutputOfCuts(Seqfasta);
  }
  fclose(Seqfasta);

  FILE *ReadSeqFile;
  ReadSeqFile=fopen(readseq_p,"w");
  OutputOfReadSeqInfo(ReadSeqFile);
  fclose(ReadSeqFile);

  //Print_OrderMatrix(parts);
  //Print_DistanceHistos(parts, 10);

  exit(0);

}  









