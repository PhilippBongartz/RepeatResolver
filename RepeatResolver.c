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

#include <gsl/gsl_cdf.h>

#include <pthread.h>

#define Max_Var_Anzahl 1500000
#define Max_Sig_Anzahl 30000

//#define PRINT



/*
Das wird die siebte Version des PrimitiveRepeatResolver
Dieses Tool soll anhand eines Multialignments Repeats auflösen.

Der Update in dieser Version besteht darin, dass ich nur Reads mit voller Coverage behandle. 
Dafür lese ich nur Reads ein, die die volle Coverage in einem bestimmten Bereich haben.
Beim Ausgeben der Unterteilung füge ich die ausgelassenen Sequenzen wieder als -1 ein.

In dieser Version versuche ich Gruppen zu bilden ohne mich auf eine Signatur zu beschränken. 

Der Algorithmus:
- Alle paarweisen CumHypGeo-Wahrscheinlichkeiten von Basengruppen werden berechnet.
- Die signifikanten Basengruppen werden mit einer Clique korrelierender Gruppen verfeinert.
- Die besten verfeinerten Gruppen werden verwendet um die Sequenzen zu unterteilen.
- Das wird rekursiv auf den Unterteilungen angewandt.
- Der Rest wird normal geclustert.

Hier gehe ich immer durch die Groups. Die mit hohen Korrs kriegen eine Clique und eine C_Group.
Außerdem einen Cutoff, eine Coverage, eine Size und eine C_Core.

*/

unsigned long *Groups[Max_Var_Anzahl*5];
unsigned long *C_Groups[Max_Var_Anzahl*5];
unsigned long *C_Coverage[Max_Var_Anzahl*5];
unsigned long *C_Core[Max_Var_Anzahl*5];
int Groupsizearray[Max_Var_Anzahl*5];
double Drop_Off[Max_Var_Anzahl*5];
int Sizes[Max_Var_Anzahl*5];
int Cutoffs[Max_Var_Anzahl*5];
int *Cliques[Max_Var_Anzahl*5];
unsigned long *LocalCoverage[Max_Var_Anzahl];
int *UnderCut;
char **Signatures;
int sc; //Diese Variable wird einmal berechnet (SigAnzahl/64 +1) und dann für das GroupHandling verwendet:
int *sc_p;  // Kann ich damit sc verändern?


void *Guarded_Malloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL)
    { fprintf(stderr,"\nError in RepeatResolver\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}

void *Guarded_Realloc(void *p, size_t size)
{ p = realloc(p,size);
  if (p == NULL)
    { fprintf(stderr,"\nError in RepeatResolver\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}


/*******************************Group-Handling**************************************/


// BitCounter geklaut von Wikipedia: Die const uint kommen in .h
const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
const uint64_t hff = 0xffffffffffffffff; //binary: all ones
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
 
//This uses fewer arithmetic operations than any other known  
//implementation on machines with fast multiplication.
//It uses 12 arithmetic operations, one of which is a multiply.
int popcount_3(uint64_t x) 
{
  x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
  x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
  x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
  return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

int intmax(int x, int y)
{
  if(x>y)return x;
  return y;
}

int intmin(int a, int b)
{
  if(a<b)return a;
  return b;
}

double doublemax(double x, double y)
{
  if(x>y)return x;
  return y;
}

double doublemin(double a, double b)
{
  if(a<b)return a;
  return b;
}


// Berechnet im Grunde die Größe des Schnitts zweier Gruppen.
int Schnitt(unsigned long *group1, unsigned long *group2)
{  
  int zz;
  int zahl=0;
  unsigned long dings;
  for (zz=0; zz<sc; zz++)
  {
    dings=*(group1+zz) & *(group2+zz);
    if(dings>0)zahl+=popcount_3(dings);                   
  }
  return zahl;  
}

// Berechnet im Grunde die Größe des Schnitts dreier Gruppen.
int Triple_Schnitt(unsigned long *group1, unsigned long *group2, unsigned long *group3)
{  
  int zz;
  int zahl=0;
  unsigned long dings;
  for (zz=0; zz<sc; zz++)
  {
    dings=*(group1+zz) & *(group2+zz) & *(group3+zz);
    if(dings>0)zahl+=popcount_3(dings);                   
  }
  return zahl;  
}

int GrMatch(unsigned long *group1, unsigned long *group2)
{
  int zz;
  int zahl=0;
  unsigned long dings;
  for (zz=0; zz<sc; zz++) //Match der 1
  {
    dings=*(group1+zz) ^ *(group2+zz);
    if(dings>0)zahl+=popcount_3(dings);                   
  }
 
  return sc*64-zahl;    
}

unsigned long *Schnitt_Group(unsigned long *group1, unsigned long *group2)
{  
  unsigned long *schnitt_group=Guarded_Malloc(sizeof(unsigned long)*sc);
  int zz;
  for (zz=0; zz<sc; zz++)
  {
    *(schnitt_group+zz) = *(group1+zz) & *(group2+zz);                   
  }
  return schnitt_group;  
}

void Inplace_Schnitt_Group(unsigned long *schnitt_group, unsigned long *group1, unsigned long *group2)
{  
  int zz;
  for (zz=0; zz<sc; zz++)
  {
    *(schnitt_group+zz) = *(group1+zz) & *(group2+zz);                   
  }
  return;
}

int Schnitt_mit_Komplement(unsigned long *group1, unsigned long *group2)
{
  int zz;
  int zahl=0;
  unsigned long dings;
  for (zz=0; zz<sc; zz++)
  {
    dings=*(group1+zz) & ~(*(group2+zz));
    zahl+=popcount_3(dings);                   
  }
  return zahl;  
}

void GrAdd(unsigned long *group, int element)
{
  unsigned long mask=1;
  int VecInd=element/64;
  mask=mask<<(element%64);
  *(group+VecInd)= *(group+VecInd) | mask;
}

void GrDel(unsigned long *group, int element)
{
  unsigned long mask=1;
  int VecInd=element/64;
  mask=mask<<(element%64);
  *(group+VecInd)= *(group+VecInd) & ~mask;
}

void GrCopy(unsigned long *group1, unsigned long *group2)
{
  int zz;
  for (zz=0; zz<sc; zz++)
  {
    *(group2+zz)= *(group1+zz);                 
  }  
}

void GrNull(unsigned long* Group)
{
  int i;
  for(i=0;i<sc;i++)
  {
    *(Group+i)=0;
  }  
}

void GrOne(unsigned long* Group)
{
  int i;
  for(i=0;i<sc;i++)
  {
    *(Group+i)=0;
    *(Group+i)=~(*(Group+i));
  }  
}

unsigned long* GrInitialize()
{
  unsigned long *group=Guarded_Malloc(sizeof(unsigned long)*sc);
  GrNull(group);
  return group;
}

int GrElement(unsigned long *group, int element)
{
  unsigned long mask=1;
  int VecInd=element/64;
  mask=mask<<(element%64);
  if(*(group+VecInd) & mask)return 1;
  return 0;
}


int Groupsize(unsigned long *group)
{
  unsigned long *castedgroup;
  int zz;
  int zahl=0;
  castedgroup=group;
  for (zz=0; zz<sc; zz++)
    {
      zahl+=popcount_3(castedgroup[zz]);               
    }
  return zahl;  
}


/******************************************Signatures-Einlesen*****************************************/
int siglength;
int signumber;
int realsigno;
int Coverage[Max_Var_Anzahl];
int Ausgelassen[Max_Sig_Anzahl];

void Einlesen(char* MApath_p, int von, int bis)   // aka Input füllen.
{ 
  char *s;  
  int i,j,k;
  char buffer[Max_Var_Anzahl];

 
  //Statt char Signatures[Max_Sig_Anzahl][Max_Var_Anzahl];  mal mallocen
  Signatures=Guarded_Malloc(sizeof(char*)*Max_Sig_Anzahl);
  for(i=0;i<Max_Sig_Anzahl;i++)
  {
    *(Signatures+i)=Guarded_Malloc(sizeof(char)*Max_Var_Anzahl);
  }

  signumber=-1;
  realsigno=-1;
  /*  //Hier lese ich vom File ein, für Debugging:
  FILE * File;
  size_t len=Max_Var_Anzahl-2;
  File=fopen("KmereMA","r");
  if(File==NULL){exit(1);}

  while ((s = fgets(buffer, len, File)) != NULL)*/

  FILE * File;
  File=fopen(MApath_p,"r");
  if(File==NULL){printf("MA is missing.\n"); exit(1);}


  while ((s=fgets(buffer, Max_Var_Anzahl-2, File)) != NULL)
  {
    realsigno++;
    siglength = strlen(buffer);
    if (buffer[siglength-1] != '\n'){printf("Was für einen Scheiß liest Du da ein? #%c# \n",buffer[siglength-1]);exit(1);}
    siglength--;  //Länge der Zeile ist inklusive '\n'
    if(bis>siglength-1){bis=siglength-1;}

    if(buffer[von]!=' ' && buffer[bis]!=' ')  
    {
      signumber++;
      Ausgelassen[realsigno]=1;
      for(i=von;i<bis+1;i++)
      {
        if(buffer[i]=='a' || buffer[i]=='A')
        {
          Signatures[signumber][i-von]=0;
        }
        else if(buffer[i]=='c' || buffer[i]=='C')
        {
          Signatures[signumber][i-von]=1;
        }
        else if(buffer[i]=='g' || buffer[i]=='G')
        {
          Signatures[signumber][i-von]=2;
        }
        else if(buffer[i]=='t' || buffer[i]=='T')
        {
          Signatures[signumber][i-von]=3;
        }
        else if(buffer[i]=='-' || buffer[i]=='_')
        {
          Signatures[signumber][i-von]=4;
        }
        else
        {
          Signatures[signumber][i-von]=5;
        }
      }      
    }

    else
    {
      Ausgelassen[realsigno]=-1;
    }
  }


  signumber++;  //Damit es nicht der letzte Index sondern die Anzahl ist.
  realsigno++;
  printf("Of %d sequences, %d had full coverage.\n",realsigno,signumber);
  printf("Siglength was %d is now %d from %d to %d.\n",siglength,bis+1-von,von,bis);
  siglength=bis+1-von; //Nur die wurden eingelesen.
  sc=(signumber/64)+1;

  #ifdef PRINT
  printf("Erfolgreich eingelesen.\n");
  fflush(stdout);
  #endif


  //Lieber Gruppen einzeln und nicht als Column organisiert sonst wird es später lästig:
  for(i=0;i<siglength*5;i++)
  {
    Groups[i]=GrInitialize();
  }  

  //Die genaue Coverage für jede Column --> SharedCoverage für zwei Gruppen.
  for(i=0;i<siglength;i++)
  {
    LocalCoverage[i]=GrInitialize();
  }  

  #ifdef PRINT
  printf("Das Allocieren ging noch.\n");
  fflush(stdout);
  #endif

  for(i=0;i<siglength;i++)Coverage[i]=0;

  for(i=0;i<siglength;i++)
  {
    for(j=0;j<signumber;j++)
    {
      for(k=0;k<5;k++)
      {
        if(Signatures[j][i]==k)
        {
          GrAdd(Groups[i*5+k],j);
        }
      }

      if(Signatures[j][i]<5)
      {
        GrAdd(LocalCoverage[i],j);
        Coverage[i]++;
      }
    }
  }
  for(i=0;i<siglength*5;i++)Groupsizearray[i]=Groupsize(Groups[i]);

  //Freeen
  for(i=0;i<Max_Sig_Anzahl;i++)
  {
    free(Signatures[i]);
  }  
  free(Signatures);
} 


double F_beta(unsigned long *Group1, unsigned long *Group2, double beta)
{
  double schnitt=(double)Schnitt(Group1, Group2); 

  //double gr1notgr2=(double)Schnitt_mit_Komplement_und_Coverage(Group1, Group2, Cov2);
  //double gr2notgr1=(double)Schnitt_mit_Komplement_und_Coverage(Group2, Group1, Cov1);

  double gr1notgr2=(double)Schnitt_mit_Komplement(Group1, Group2);
  double gr2notgr1=(double)Schnitt_mit_Komplement(Group2, Group1);

  double Z=(1.0+beta)*(double)schnitt;
  if(Z<0.0001)return 0.0;
  Z/=( (1+beta*beta)*schnitt+(beta*beta*gr1notgr2)+gr2notgr1);

  return Z; 
}

double PositiveCumHypGeo_Log(unsigned int schnitt, unsigned int gr1, unsigned int gr2, unsigned int cov)
{ 
  double posQ=gsl_cdf_hypergeometric_Q (schnitt-1, gr2, cov-gr2, gr1); 
  posQ=-1.0*log10(posQ);
  if(isinf(posQ) || posQ>99)return 99.0;
  return posQ;
}

double PositiveSignificance(int i, int j)
{
  int schnitt=Schnitt(Groups[i], Groups[j]);   
  int cov=Schnitt(LocalCoverage[i/5],LocalCoverage[j/5]); 
  int gr1=Schnitt(Groups[i],LocalCoverage[j/5]);
  int gr2=Schnitt(Groups[j],LocalCoverage[i/5]);
  if(gr1==0 || gr2==0)return 0.0;
  //if(schnitt<gr1*gr2/cov*2 || schnitt<1)return 0.0; //Eine Bremse
  if(schnitt<1)return 0.0; //Eine Bremse
  double Z=PositiveCumHypGeo_Log(schnitt,gr1,gr2,cov);  
  if(isinf(Z) || Z>98.0)Z=98.0+F_beta(Groups[i], Groups[j], 1.0);
  return Z; 
}


double Group_PositiveSignificance(unsigned long *Group1,unsigned long *Group2,unsigned long *Cov1,unsigned long *Cov2)
{
  //printf("Drin\n");fflush(stdout);
  int schnitt=Schnitt(Group1, Group2);   
  //printf("S%d\n",schnitt );fflush(stdout);
  int cov=Schnitt(Cov1,Cov2); 
  //printf("C%d\n",cov );fflush(stdout);
  int gr1=Schnitt(Group1,Cov2);
  //printf("g1 %d\n",gr1 );fflush(stdout);
  int gr2=Schnitt(Group2,Cov1);
  if(gr1==0 || gr2==0)return 0.0;
  //printf("g2 %d\n",gr2 );fflush(stdout);
  double Z=PositiveCumHypGeo_Log(schnitt,gr1,gr2,cov);  
  //printf("%f \n\n",Z );
  if(isinf(Z) || Z>98.0)Z=97.90+F_beta(Group1, Group2, 1.0);
  return Z; 
}

double CumHypGeo_Log(unsigned int schnitt, unsigned int gr1, unsigned int gr2, unsigned int cov)
{ 
  double posP=gsl_cdf_hypergeometric_P (schnitt, gr2, cov-gr2, gr1);  
  double posQ=gsl_cdf_hypergeometric_Q (schnitt-1, gr2, cov-gr2, gr1); 

  if(posP<posQ || schnitt==0)
  {
    posP=-1.0*log10(posP);
    if(isinf(posP) || posP>99)return 99.0;
    return posP;
  }
  posQ=-1.0*log10(posQ);
  if(isinf(posQ) || posQ>99)return 99.0;
  return posQ;
}

double Relative_Group_Significance(unsigned long *Group1,unsigned long *Group2,unsigned long *Cov)
{
  //printf("Drin\n");fflush(stdout);
  int schnitt=Triple_Schnitt(Group1, Group2, Cov);   
  //printf("S%d\n",schnitt );fflush(stdout);
  int cov=Groupsize(Cov); 
  //printf("C%d\n",cov );fflush(stdout);
  int gr1=Schnitt(Group1,Cov);
  //printf("g1 %d\n",gr1 );fflush(stdout);
  int gr2=Schnitt(Group2,Cov);
  if(gr1==0 || gr2==0)return 0.0;
  //printf("g2 %d\n",gr2 );fflush(stdout);
  double Z=CumHypGeo_Log(schnitt,gr1,gr2,cov);  
  //printf("%f \n\n",Z );
  if(isinf(Z) || Z>99.0)Z=99.0;
  return Z; 
}

double Erwartete_Anzahl_von_Sigs_mit_mehr_als_c_positives(int n,int c, double p,int v)
{
  return gsl_cdf_binomial_Q(c,p,v)*n; //k,p,n
}


int BestCutoff(int n, int nn, int v, double p, double pp)
{
  int c,bestc;
  bestc=0;
  double score,bestscore;
  bestscore=0;
  for(c=0;c<v;c++)
  {
    score=Erwartete_Anzahl_von_Sigs_mit_mehr_als_c_positives(n,c,p,v);
    score/=doublemax(Erwartete_Anzahl_von_Sigs_mit_mehr_als_c_positives(nn,c,pp,v),1.0);
    //printf("%d/%d -> %f=%f/%f\n",c,v,score,Erwartete_Anzahl_von_Sigs_mit_mehr_als_c_positives(n,c,p,v),Erwartete_Anzahl_von_Sigs_mit_mehr_als_c_positives(nn,c,pp,v));
    if(score>bestscore)
    {
      bestscore=score;
      bestc=c;
    }
  }
  return bestc;
}  

void MaxCorrsRausschreiben(double *MaxCorrs, char *outputfile)
{
  FILE *datei;
  datei=fopen(outputfile,"w");

  if(NULL == datei)
  {
    printf("DateiVerbratei!\n");
    exit(1);
  }
  int i;
  for(i=0;i<siglength*5;i++)
  {
    fprintf(datei,"%f\n", MaxCorrs[i]);
  }
  fclose(datei);
}

void Unterteilung_Rausschreiben(int *Unterteilung, char *outputfile)
{
  FILE *datei;
  datei=fopen(outputfile,"w");

  if(NULL == datei)
  {
    printf("DateiVerbratei!\n");
    exit(1);
  }
  int i;
  for(i=0;i<realsigno;i++)  //Nicht signumber, das kann nur eine Teilmenge sein. 
  {
    if(i!=0){fprintf(datei,"\n");}
    fprintf(datei,"%d", Unterteilung[i]);
  }
  fclose(datei);
}

int *UnterteilungEinlesen(char *inputfile)
{
  char *s;
  FILE * File;
  size_t len=100;
  File=fopen(inputfile,"r");
  if(File==NULL){fclose(File);return NULL;}
  char buffer[101];

  int *Unterteilung=Guarded_Malloc(sizeof(int)*signumber);
  int i=0;

  while ((s = fgets(buffer, len, File)) != NULL)
  {
    sscanf(&buffer[0],"%d",(Unterteilung+i));
    i++;
  }
  return Unterteilung;
}



double *MaxCorrsEinlesen(char *inputfile, int von, int bis)
{
  char *s;
  FILE * File=NULL;
  size_t len=100;
  File=fopen(inputfile,"r");

  #ifdef PRINT
  printf("%p\n",File );
  fflush(stdout);
  #endif

  if(File==NULL){return NULL;}

  char buffer[101];

  double *MaxCorrs=Guarded_Malloc(sizeof(double)*siglength*5);
  int i=0;
  int j=0;

  while ((s = fgets(buffer, len, File)) != NULL)
  {
    if(i/5>=von && i/5<=bis)
    {
      sscanf(&buffer[0],"%lf",(MaxCorrs+j));
      //printf("%f\n", *(MaxCorrs+j));
      j++;
    }
    i++;
  }

  #ifdef PRINT  
  printf("%d MaxCorrs eingelesen %d - %d\n",j,von,bis );
  #endif

  fclose(File);
  return MaxCorrs;
}


double *AllMaxCorrsRechner(int anfang, int ende, int mincov, int maxgroup, double cutoff)
{

  int i,ii,j,jj,k,kk;
  int *count=Guarded_Malloc(sizeof(int)*siglength*5);
  for(i=0;i<siglength*5;i++)count[i]=0;
  //int *Groupsizearray=Guarded_Malloc(sizeof(int)*siglength*5);
  //for(i=0;i<siglength*5;i++)Groupsizearray[i]=Groupsize(Groups[i]); 

  double *MaxCorrs=Guarded_Malloc(sizeof(double)*siglength*5);
  for(i=0;i<siglength*5;i++){MaxCorrs[i]=0.0;}

  double Z;
  //Erstmal suche ich die beiden am stärksten korrelierenden relevanten Vars raus:
  for(ii=anfang;ii<ende;ii++)
  {
    //printf("%d\n",ii);fflush(stdout);
    for(k=0;k<5;k++)
    {
      i=ii*5+k; //Die relevante Var in Column ii
      if(Groupsizearray[i]>mincov/4 && Groupsizearray[i]<maxgroup)
      {
        for(jj=ii+20;jj<ende;jj++)  
        {
          //Wenn der Abstand zu groß wird gibt es keine gemeinsame coverage mehr, dann wird abgebrochen:
          if(Schnitt(LocalCoverage[ii],LocalCoverage[jj])<mincov)
          {
            jj=ende;
          }

          else
          {
            for(kk=0;kk<5;kk++)
            {
              j=jj*5+kk; //Die relevante Var in Column jj
              if(Groupsizearray[j]>mincov/4 && Groupsizearray[j]<maxgroup)
              {
                Z=PositiveSignificance(i,j);
                if(Z>MaxCorrs[i]){MaxCorrs[i]=Z;}
                if(Z>MaxCorrs[j]){MaxCorrs[j]=Z;}
                if(Z>cutoff){count[i]++;count[j]++;}
              }
            }
          }
        } 
      }
    }
  } 
  //Hier kicke ich die raus, die zu selten korrelieren:
  for(i=0;i<siglength*5;i++)
  {
    if(count[i]<5)MaxCorrs[i]=0.0;
  }
  return MaxCorrs;
}  

double *TestCorrsRechner(int anfang, int ende, int mincov, int maxgroup, double cutoff)
{

  int i,j,jj,ii;
  //int *Groupsizearray=Guarded_Malloc(sizeof(int)*siglength*5);
  //for(i=0;i<siglength*5;i++)Groupsizearray[i]=Groupsize(Groups[i]); 

  double *MaxCorrs=Guarded_Malloc(sizeof(double)*siglength*5);
  for(i=0;i<siglength*5;i++){MaxCorrs[i]=0.0;}

  int histo[50];
  for(i=0;i<50;i++)histo[i]=0;

  int count;
  double Z;
  //Erstmal suche ich die beiden am stärksten korrelierenden relevanten Vars raus:
  for(i=anfang;i<ende*5;i++)
  {
    count=0;
    if(Groupsizearray[i]>mincov/4 && Groupsizearray[i]<maxgroup)
    {
      for(j=anfang;j<ende*5;j++)  
      {
        if(abs(i-j)>100)
        {
          if(Groupsizearray[j]>mincov/4 && Groupsizearray[j]<maxgroup)
          {
            Z=PositiveSignificance(i,j);
            if(Z>MaxCorrs[i]){MaxCorrs[i]=Z;}
            if(Z>MaxCorrs[j]){MaxCorrs[j]=Z;}
            if(Z>cutoff){count++;jj=j;}
          }
        }
      } 
    }
    histo[intmin(49,count)]++;
    #ifdef PRINT
    for(ii=0;ii<50;ii++)printf("%d ",histo[ii]);
    printf("\n");
    if(count==1 && MaxCorrs[i]>50)printf("i:%d,j=%d,gr1=%d,gr2=%d,Z=%f,\n", i,jj,Groupsizearray[i],Groupsizearray[jj],MaxCorrs[i]);
    #endif
  } 
  return MaxCorrs;
} 

void CheckAll(int anfang, int ende, int mincov, int maxgroup, double cutoff, int i, double *MaxCorrs, int *Used, int *Groupsizearray)
{
  if(Used[i])return;
  Used[i]=1;
  //printf("%d\n",i );

  int jj,j,kk;
  double Z;
  if(Groupsizearray[i]>mincov/4 && Groupsizearray[i]<maxgroup)
  {
    for(jj=anfang;jj<ende;jj++)  
    {
      for(kk=0;kk<5;kk++)
      {
        j=jj*5+kk; //Die relevante Var in Column jj
        if(Groupsizearray[j]>mincov/4 && Groupsizearray[j]<maxgroup && Used[j]==0 && abs(j-i)>100)
        {
          Z=PositiveSignificance(i,j);
          if(Z>MaxCorrs[i]){MaxCorrs[i]=Z;}
          if(Z>MaxCorrs[j])
          {
            MaxCorrs[j]=Z;
            if(Z>cutoff && Used[j]==0)
            {
              CheckAll(anfang, ende, mincov, maxgroup, cutoff, j, MaxCorrs, Used, Groupsizearray);
            }
          }
        }
      }
    } 
  }  
}

double *Snowball_AllMaxCorrsRechner(int anfang, int ende, int mincov, int maxgroup, double cutoff,int fraction)
{

  int i,ii,k;
  //int *Groupsizearray=Guarded_Malloc(sizeof(int)*siglength*5);
  //for(i=0;i<siglength*5;i++)Groupsizearray[i]=Groupsize(Groups[i]); 

  double *MaxCorrs=Guarded_Malloc(sizeof(double)*siglength*5);
  for(i=0;i<siglength*5;i++){MaxCorrs[i]=0.0;}

  int *Used=Guarded_Malloc(sizeof(int)*siglength*5);
  for(i=0;i<siglength*5;i++)Used[i]=0;

  //Erstmal suche ich die beiden am stärksten korrelierenden relevanten Vars raus:
  for(ii=anfang;ii<ende;ii++)
  {
    //printf("%d\n",ii);fflush(stdout);
    for(k=0;k<5;k++)
    {
      i=ii*5+k; //Die relevante Var in Column ii
      if(Groupsizearray[i]>signumber/fraction && Groupsizearray[i]<signumber-signumber/fraction && Used[i]==0)  //erster Test: Nur die großen Groups.
      {
        CheckAll(anfang, ende, mincov, maxgroup, cutoff, i, MaxCorrs, Used, Groupsizearray);
      }
    }
  } 
  return MaxCorrs;
}  

double *HilfsMaxCorr[64];
void *HilfsMaxCorrsRechner(void *x)
{
  time_t start_time=time(NULL);

  int * args;
  args=((int *) x);

  int anfang=(*(args+0));
  int ende=(*(args+1)); 
  int mincov=(*(args+2));
  int maxgroup=(*(args+3)); 
  double cutoff=(*(args+4)); 
  int NTHREADS=(*(args+5)); 
  int thread=(*(args+6)); 

  double *MaxCorrs;

  int i,ii,j,jj,k,kk;
  int *count=Guarded_Malloc(sizeof(int)*siglength*5);
  for(i=0;i<siglength*5;i++)count[i]=0;

  MaxCorrs=Guarded_Malloc(sizeof(double)*siglength*5);
  HilfsMaxCorr[thread]=MaxCorrs;
  for(i=0;i<siglength*5;i++){MaxCorrs[i]=0.0;}

  double Z;
  int percentage=1;
  time_t current_time;
  //Erstmal suche ich die beiden am stärksten korrelierenden relevanten Vars raus:
  for(ii=anfang;ii<ende;ii++)
  {
    if(thread==0)
    {
      //printf("%d ",ii );
      if( (percentage*(ende-anfang))/100 == ii )
      {
        current_time=time(NULL);
        #ifdef PRINT
        //printf("%d %% in %ld min -> %ld min.\n",percentage,(current_time-start_time)/60,(((current_time-start_time)*(100-percentage))/(percentage*60))  );  
        printf("%d %% in %ld min of %ld min\n",percentage,(current_time-start_time)/60, (10000*(current_time-start_time)/ ( (200-percentage)*percentage ))/60);
        #endif
        percentage++;
      }
      fflush(stdout); 
    }     


    if(ii%NTHREADS==thread)
    {
      for(k=0;k<5;k++)
      {
        i=ii*5+k; //Die relevante Var in Column ii
        if(Groupsizearray[i]>mincov/4 && Groupsizearray[i]<maxgroup)
        {
          for(jj=ii+20;jj<ende;jj++)  
          {
            //Wenn der Abstand zu groß wird gibt es keine gemeinsame coverage mehr, dann wird abgebrochen:
            if(Schnitt(LocalCoverage[ii],LocalCoverage[jj])<mincov)
            {
              jj=ende;
            }

            else
            {
              for(kk=0;kk<5;kk++)
              {
                j=jj*5+kk; //Die relevante Var in Column jj
                if(Groupsizearray[j]>mincov/4 && Groupsizearray[j]<maxgroup)
                {
                  Z=PositiveSignificance(i,j);
                  if(Z>MaxCorrs[i]){MaxCorrs[i]=Z;}
                  if(Z>MaxCorrs[j]){MaxCorrs[j]=Z;}
                  if(Z>cutoff){count[i]++;count[j]++;}
                }
              }
            }
          } 
        }
      }
    }
  } 
  #ifdef PRINT
  printf("thread %d done\n",thread );fflush(stdout);
  #endif
  return NULL;
}  

double *Parallel_AllMaxCorrsRechner(int NTHREADS, int anfang, int ende, int mincov, int maxgroup, double cutoff)
{

  time_t start_time=time(NULL);

  pthread_t threads[NTHREADS];
  int *thread_args[NTHREADS];

  int rc,i;

  /* spawn the threads */
  for (i=0; i<NTHREADS; ++i)
  {
    //(int anfang, int ende, int mincov, int maxgroup, double cutoff, int NTHREADS, int thread)
    thread_args[i] = Guarded_Malloc(sizeof(int)*7);
    thread_args[i][0]=anfang;
    thread_args[i][1]=ende;
    thread_args[i][2]=mincov;
    thread_args[i][3]=maxgroup;
    thread_args[i][4]=cutoff;
    thread_args[i][5]=NTHREADS;
    thread_args[i][6]=i;

    #ifdef PRINT
    printf("spawning thread %d\n", i);
    #endif
    rc = pthread_create(&threads[i], NULL, HilfsMaxCorrsRechner, (void *) thread_args[i]);
  }

  #ifdef PRINT
  printf("threads done\n");fflush(stdout);
  #endif

  /* wait for threads to finish */
  for (i=0; i<NTHREADS; ++i) {
    rc = pthread_join(threads[i], NULL);
  }

  #ifdef PRINT
  printf("pthread_join\n");fflush(stdout);
  #endif

  // MaxCorrs vereinigen:
  for(i=0;i<siglength*5;i++)
  {
    HilfsMaxCorr[0][i]=doublemax(HilfsMaxCorr[0][i],doublemax(HilfsMaxCorr[1][i],HilfsMaxCorr[2][i]));
  }

  #ifdef PRINT
  printf("MaxCorrs vereinigt\n");fflush(stdout);
  #endif

  //clock_t end = clock();
  //double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  time_t time_spent=time(NULL)-start_time;
  #ifdef PRINT
  printf("%ld sec.\n", time_spent);  
  #endif

  for (i=1; i<NTHREADS; ++i)free(HilfsMaxCorr[i]);

  return HilfsMaxCorr[0];
}











// Berechnet eine Group für eine Clique und einen gegebenen Cutoff
unsigned long *CliqueGroup(int *Clique, int c) //, int sig)
{
  unsigned long *Group=Guarded_Malloc(sizeof(unsigned long)*sc);
  GrNull(Group);

  int i,ii,jj,j;
  for(jj=0;jj<100;jj++)  //Die Größe der Clique feststellen.
  {
    if(Clique[jj]<0)
    {
      j=jj;
      jj=100;
    }
  }

  for(i=0;i<signumber;i++)
  {
    ii=0; //Der score
    for(jj=0;jj<j;jj++) //Durch die vars
    {
      //if(Signatures[i][Clique[jj]/5]==Signatures[sig][Clique[jj]/5])
      if(GrElement(Groups[Clique[jj]],i))
      {
        ii++;
      }
    }
    if(ii>c)
    {
      GrAdd(Group,i);
    }
  }
  return Group;
}

// Berechnet eine Group für eine Clique und einen gegebenen GruppenGrößenCutoff! 
unsigned long *CoreGroup(int *Clique, int c) //, int sig)
{
  unsigned long *Group=Guarded_Malloc(sizeof(unsigned long)*sc);
  GrNull(Group);

  int i,ii,jj,j;
  for(jj=0;jj<100;jj++)  //Die Größe der Clique feststellen.
  {
    if(Clique[jj]<0)
    {
      j=jj;
      jj=100;
    }
  }
  int histo[100];
  for(i=0;i<j;i++)histo[i]=0;
  for(i=0;i<signumber;i++)
  {
    ii=0; //Der score
    for(jj=0;jj<j;jj++) //Durch die vars
    {
      //if(Signatures[i][Clique[jj]/5]==Signatures[sig][Clique[jj]/5])
      if(GrElement(Groups[Clique[jj]],i))
      {
        histo[ii]++;
        ii++;
      }
    }
  }
  //Cutoff feststellen:
  i=0;
  while(histo[i]>c)i++;
  c=i;
  //Elemente hinzufügen:
  for(i=0;i<signumber;i++)
  {
    ii=0; //Der score
    for(jj=0;jj<j;jj++) //Durch die vars
    {
      //if(Signatures[i][Clique[jj]/5]==Signatures[sig][Clique[jj]/5])
      if(GrElement(Groups[Clique[jj]],i))
      {
        ii++;
      }
    }
    if(ii>c)
    {
      GrAdd(Group,i);
    }
  }
  return Group;
}

unsigned long *CliqueCoverage(int *Clique, int c) //, int sig)
{
  unsigned long *Coverage=Guarded_Malloc(sizeof(unsigned long)*sc);
  GrNull(Coverage);

  int i,ii,jj,j;
  for(jj=0;jj<100;jj++)  //Die Größe der Clique feststellen.
  {
    if(Clique[jj]<0)
    {
      j=jj;
      jj=100;
    }
  }

  for(i=0;i<signumber;i++)
  {
    ii=0; //Der score
    for(jj=0;jj<j;jj++) //Durch die vars
    {
      //if(Signatures[i][Clique[jj]/5]<5)
      if(GrElement(LocalCoverage[Clique[jj]/5],i))  //Irgendwann alles ohne Signatures.
      {
        ii++;
      }
    }
    if(ii>c)
    {
      GrAdd(Coverage,i);
    }
  }
  return Coverage;
}

void GroupPrecision(unsigned long *group)
{
  int i,j;
  int maj=0;
  int min=0;
  int drin=0;
  int drau=0;
  for(i=0;i<signumber/30;i++)
  {
    drin=0;
    drau=0;
    for(j=0;j<30;j++)
    {
      if(GrElement(group,i*30+j))
      {drin++;}
      else
      {
        drau++;
      }
    }
    if(drin>drau)
    {
      maj+=drin;
      min+=drau;
    }
    else
    {
      maj+=drau;
      min+=drin;
    }
  }

  printf("Group Precision %d / %d\n",maj,min);
}

void GroupPrinter(unsigned long *group)
{
  int i,j;
  int drin=0;
  int drau=0;
  for(i=0;i<signumber/30;i++)
  {
    drin=0;
    drau=0;
    for(j=0;j<30;j++)
    {
      if(GrElement(group,i*30+j))
      {drin++;}
      else
      {
        drau++;
      }
    }
    printf("%d ",drin );
  }
  printf("\n");
}

void TheBestUpdater(int *Clique, double *Best_Corrs, int maxclique, int i, double Z)
{
  if(Best_Corrs[maxclique-1]>=Z){return;} //Keine Verbesserung möglich.

  int ii,j;
  ii=maxclique-1; //Letzter Eintrag
  while(Best_Corrs[ii]<Z && ii>0) //Neuer Eintrag ist besser
  {
    ii--;
  }
  ii++; //Da gehört der neue Eintrag hin alle anderen werden nach hinten verschoben.

  for(j=maxclique-1;j>ii;j--)
  {
    Best_Corrs[j]=Best_Corrs[j-1];
    Clique[j]=Clique[j-1];
  }

  Best_Corrs[ii]=Z;
  Clique[ii]=i;
}

// Hier wird eine Clique greedy nur auf der ersten Gruppe aufgebaut.
int *Cliquer(int anfang, int ende, int mincov, int maxclique, double greedy, int a)
{

  #ifdef PRINT
  printf("\nCliquer\n");
  #endif

  int *Clique=Guarded_Malloc(sizeof(int)*(maxclique+1));
  double *Best_Corrs=Guarded_Malloc(sizeof(double)*(maxclique+1));
  int schnitt; //, cov, gr1,gr2;
  int ii,i,j,k;//,jj;

  //int gr=Groupsize(Groups[a]);

  double Z;
  //double maxcorr=greedy+1.0;

  Clique[0]=a;
  for(i=0;i<maxclique;i++)Best_Corrs[i]=0.0;

  j=1; //Wird jetzt als Cliquengröße verwendet.

  //printf("%f %f %d %d\n",maxcorr,greedy,j,maxclique );

  //Hier fülle ich die Clique mit dem BestUpdater:
  for(ii=anfang;ii<ende;ii++) //Für alle columns
  {
    for(k=0;k<5;k++)
    {
      i=ii*5+k; //Die relevante Var/Group in Column ii
      if(i!=Clique[0])
      {
        Z=0.0;
        schnitt=Schnitt(Groups[i], Groups[Clique[0]]);   
        fflush(stdout);
        if(schnitt>mincov/4)// && ((double)(abs(Groupsize(Groups[i])-gr)))/((double)gr)<0.5 )
        {
          Z=Group_PositiveSignificance(Groups[i],Groups[Clique[0]],LocalCoverage[ii],LocalCoverage[Clique[0]/5]); 
          if(Z>greedy)
          {
            TheBestUpdater(Clique, Best_Corrs, maxclique, i, Z);
            //for(i=0;i<maxclique;i++)printf("%f ",Best_Corrs[i]);
            //printf("\n");
          }
        }  
      }
    }
  }

  Best_Corrs[0]=100.0;
  Clique[maxclique]=-1;
  j=maxclique-1;
  while(Best_Corrs[j]<greedy || Clique[j]==Clique[j-1]){Clique[j]=-1;j--;}

  #ifdef PRINT
  for(i=0;i<maxclique;i++)printf("%f ",Best_Corrs[i]);
  printf("\n");
  fflush(stdout);
  #endif

  return Clique;
}

// Hier wird eine Clique greedy auf einem Paar aufgebaut, aber die corrs werden mit der Clique berechnet.
// Nicht mit den Elementen der Clique.
int *Cliquer2(int anfang, int ende, int mincov, int maxclique, double greedy, int a)
{

  #ifdef PRINT
  printf("\nCliquer2\n");
  fflush(stdout);
  #endif

  int *Clique=Guarded_Malloc(sizeof(int)*100);
  int schnitt, cov, gr1,gr2;
  int ii,i,jj,j,k,c;
  double maxcorr=100.0;
  double currentcorr,Z;

  int n=30;  //Das muss durch eine Schätzung ersetzt werden
  int nn=signumber;
  double p=0.70;
  double pp=0.05;

  for(i=0;i<100;i++){Clique[i]=-1;} //Damit man das Ende erkennt.
  Clique[0]=a;

  long unsigned *C_Group=CliqueGroup(Clique,0);
  long unsigned *C_Coverage=CliqueCoverage(Clique,0);

  j=1; //Wird jetzt als Cliquengröße verwendet.

  //Jetzt wird die Clique greedy gefüllt:
  while(maxcorr>greedy && j<maxclique+1)
  {
    maxcorr=0;
    for(ii=anfang;ii<ende;ii++) //Für alle columns
    {
      for(k=0;k<5;k++)
      {
        i=ii*5+k; //Die relevante Var/Group in Column ii
        currentcorr=0.0;

        for(jj=0;jj<j;jj++)  //Wird gecheckt ob man die Var schon verwendet:
        {
          if(Clique[jj]==i)currentcorr-=1000.0;
        }

        schnitt=Schnitt(Groups[i], C_Group);   
        cov=Schnitt(LocalCoverage[ii],C_Coverage); 
        gr1=Schnitt(Groups[i],C_Coverage);
        gr2=Schnitt(C_Group,LocalCoverage[ii]);
        Z=0.0;
        //printf("%d %d %d %d %f\n",schnitt,cov,gr1,gr2,Z );
        fflush(stdout);
        if(schnitt>mincov/4 && gr1<500)Z=Group_PositiveSignificance(Groups[i],C_Group,LocalCoverage[ii],C_Coverage); //PositiveCumHypGeo_Log(schnitt,gr1,gr2,cov);
        //printf("%d %d %d %d %f\n",schnitt,cov,gr1,gr2,Z );
        currentcorr+=Z;

        if(currentcorr>maxcorr)
        {
          Clique[j]=i;  //vorläufiges Cliquenmitglied wird evtl überschrieben
          maxcorr=currentcorr;
        }
      }
    }
    #ifdef PRINT
    printf("%d - %d = %f\n",j,Clique[j],maxcorr );
    #endif

    j++;
    free(C_Group);
    free(C_Coverage);
    C_Group=CliqueGroup(Clique,c);
    C_Coverage=CliqueCoverage(Clique,c);  
    c=BestCutoff(n,nn,j,p,pp);  

  }//end while
  j--; //Letztes ist kleiner als greedy
  Clique[j]=-1;  //Als stopsign.

  return Clique;
}

// Hier wird eine Clique greedy auf einem Paar aufgebaut:
int *Cliquer3(int anfang, int ende, int mincov, int maxclique, double greedy, int a)
{

  #ifdef PRINT
  printf("\nCliquer3\n");
  #endif

  int *Clique=Guarded_Malloc(sizeof(int)*100);
  int schnitt,cov,gr1,gr2;
  int ii,i,jj,j,k;
  double maxcorr=100.0;
  double currentcorr,Z;

  Clique[0]=a;
  #ifdef PRINT
  printf("Groups %d and %d corr %f\n",Clique[0],Clique[1],maxcorr);
  #endif

  j=1; //Wird jetzt als Cliquengröße verwendet.

  //Jetzt wird die Clique greedy gefüllt:
  while(maxcorr>greedy && j<maxclique)
  {
    maxcorr=0;
    for(ii=anfang;ii<ende;ii++) //Für alle columns
    {
      for(k=0;k<5;k++)
      {
        i=ii*5+k; //Die relevante Var/Group in Column ii
        currentcorr=0.0;
        for(jj=0;jj<j;jj++)  //Wird die durchschnittliche Corr mit den Cliquenmitgliedern berechnet
        {
          if(Clique[jj]==i)currentcorr-=1000.0;
          else
          {
            //printf("%d %d %d %d \n",i,Clique[jj],jj,siglength*5 );
            schnitt=Schnitt(Groups[i], Groups[Clique[jj]]);   
            cov=Schnitt(LocalCoverage[ii],LocalCoverage[Clique[jj]/5]); 
            gr1=Schnitt(Groups[i],LocalCoverage[Clique[jj]/5]);
            gr2=Schnitt(Groups[Clique[jj]],LocalCoverage[ii]);
            Z=0.0;
            //printf("%d %d %d %d %f\n",schnitt,cov,gr1,gr2,Z );
            fflush(stdout);
            if(schnitt>mincov/4)Z=PositiveCumHypGeo_Log(schnitt,gr1,gr2,cov);
            //printf("%d %d %d %d %f\n",schnitt,cov,gr1,gr2,Z );
            currentcorr+=Z;
          }
        }
        currentcorr/=(double)j;
        if(currentcorr>maxcorr)
        {
          Clique[j]=i;
          maxcorr=currentcorr;
        }
      }
    }
    #ifdef PRINT
    printf("%d - %d = %f\n",j,Clique[j],maxcorr );
    #endif
    j++;

  }//end while
  j--; //Letztes ist kleiner als gr
  Clique[j]=-1;  //Als stopsign.

  return Clique;
}

// Hier eine Funktion, die den Cutoff sucht, der die Korrelation der C_Group mit der Group maximiert
int KorrMaxCutoff(int c_i, int c)
{

  unsigned long *Cutoff_Groups[100];
  int i,jj,ii;
  for(i=0;i<Sizes[c_i];i++)
  {
    Cutoff_Groups[i]=GrInitialize();
  }

  for(i=0;i<signumber;i++)
  {
    ii=0; //Der score
    for(jj=0;jj<Sizes[c_i];jj++) //Durch die vars
    {
      if(GrElement(Groups[Cliques[c_i][jj]],i))
      {
        GrAdd(Cutoff_Groups[ii],i);  //In mehr als ii Groups
        ii++;
      }
    }
  }  

  float possig;
  float maxps=0.0;
  int max_i=0;

  for(i=c;i<Sizes[c_i];i++)  //Größer als der vorherige Cutoff, der die Anzahl von fp auf unter 1 drückt.
  {
    possig=Group_PositiveSignificance(Cutoff_Groups[i],Groups[Cliques[c_i][0]],LocalCoverage[Cliques[c_i][0]/5],LocalCoverage[Cliques[c_i][0]/5]);
    #ifdef PRINT
    printf("%f ", possig);  
    #endif
    //GroupPrecision(Cutoff_Groups[i]);
    if(possig>maxps)
    {
      maxps=possig;
      max_i=i;
    }
  }  

  #ifdef PRINT
  printf("\n");

  for(i=0;i<Sizes[c_i];i++)
  {
    printf("%d ", Groupsize(Cutoff_Groups[i]));
  }  
  printf("\n");  

  for(i=0;i<Sizes[c_i];i++)
  {
    printf("%d ", Schnitt(Cutoff_Groups[i],Groups[Cliques[c_i][0]]));
  }  
  printf("\n");   
  printf("Grsize: %d\n",Groupsize(Groups[Cliques[c_i][0]]) ); 
  #endif

  for(i=0;i<Sizes[c_i];i++)
  {
    free(Cutoff_Groups[i]);
  }

  return max_i;
}

// Der Cutoff, der den Dropoff der C_Groupsizes minimiert:
int Dropoff_Cutoff(int c_i, int c)
{

  unsigned long *Cutoff_Groups[100];
  int i,jj,ii;

  // Die Groups für verschiedene Cutoffs:
  for(i=0;i<Sizes[c_i];i++)
  {
    Cutoff_Groups[i]=GrInitialize();
  }
  for(i=0;i<signumber;i++)
  {
    ii=0; //Der score
    for(jj=0;jj<Sizes[c_i];jj++) //Durch die vars
    {
      if(GrElement(Groups[Cliques[c_i][jj]],i))
      {
        GrAdd(Cutoff_Groups[ii],i);  //In mehr als ii Groups
        ii++;
      }
    }
  }  

  // Minimalen Dropoff:
  double sizes[100];
  for(i=0;i<Sizes[c_i];i++)sizes[i]=(double)Groupsize(Cutoff_Groups[i]);
  int drop_c=intmax(1,c);
  c=drop_c;
  double min_drop=1000000.0; 
  double drop;
  #ifdef PRINT
  printf("Drops %d: ",c);
  #endif
  for(i=c;i<Sizes[c_i]-1;i++)
  {
    if(doublemin((double)signumber-sizes[i],sizes[i])>0)
    {
      drop=((sizes[i-1]-sizes[i+1])/doublemin((double)signumber-sizes[i],sizes[i]));
      #ifdef PRINT
      printf("%f ",drop);
      #endif
      if(drop<min_drop)
      {
        min_drop=drop;
        drop_c=i;
      }
    }
  }
  #ifdef PRINT
  printf("\n");
  #endif

  // Den Wert eintragen.
  Drop_Off[c_i]=min_drop;

  for(i=0;i<Sizes[c_i];i++)
  {
    free(Cutoff_Groups[i]);
  }

  return drop_c;
}

// Der Cutoff, der den Dropoff der C_Groupsizes minimiert und zwar auch über Sizes:
int Sizes_Dropoff_Cutoff(int c_i, int c, int minsize)
{
  int i,k,s,count;
  int **Results=Guarded_Malloc(sizeof(int*)*100);
  for(i=0;i<100;i++)Results[i]=Guarded_Malloc(sizeof(int)*100);
  for(i=0;i<100;i++)
  {
    for(k=0;k<100;k++){Results[i][k]=0;}
  }

  for(i=0;i<signumber;i++)
  {
    count=0;
    for(k=0;k<Sizes[c_i];k++)
    {
      if(GrElement(Groups[Cliques[c_i][k]],i)){count++;}

      for(s=0;s<count;s++)Results[k+1][s]++;    //For Cliquesize k+1 and cutoff < count ist sig in consensus
    }
  }

  // Minimalen Dropoff:
  int drop_c=intmax(1,c);
  c=drop_c;

  double min_drop=1000000.0; 
  int min_size=0;
  double drop;

  //printf("Drops %d: ",c);
  for(k=minsize;k<Sizes[c_i]+1;k++)  //For different sizes k
  {
    //printf("\nSize%d: ",k);
    for(i=c;i<k;i++) // For cutoffs kleiner als size k
    {
      if(doublemin((double)signumber-Results[k][i],Results[k][i])>0) //If there aren't zero
      {
        drop=((Results[k][i-1]-Results[k][i+1])/doublemin((double)signumber-Results[k][i],Results[k][i]));
        //printf("%f ",drop);
        if(drop<=min_drop)
        {
          min_drop=drop;
          drop_c=i;
          min_size=k;
        }
      }
    }
  }

  //printf("\n");

  // Den Wert eintragen.
  Drop_Off[c_i]=min_drop;
  Sizes[c_i]=min_size;

  for(i=0;i<100;i++)free(Results[i]);
  free(Results);

  return drop_c;
}


//Für gr<100 p=0.77, gr<1000 p=83, 1500<gr<2000 p=0.90, gr >2500 p=93.0
//Das muss ich noch durch anständig intrapolierte Werte ersetzen. 
int Core_Cutoff(int gr)
{
  double p;
  if(gr<100)p=0.77;
  else if(gr<1000)p=0.83;
  else if(gr<2000)p=0.90;
  else p=0.93;
  return (int)(p*gr);
}


//Bei wenigen Vars ersetze ich Groups durch Schnitte um präzisere Auflösung zu erreichen.
// Inplace_Schnitt_Group(unsigned long *schnitt_group, unsigned long *group1, unsigned long *group2)
// Entweder alle durch Schnitte ersetzen oder Schnitte hinzufügen um normale Groups zu refinen.
// Das Problem ist die LocalCoverage ... 
void Cut_and_Paste(double *MaxCorrs, double cutoff, int mincov)
{
  int i,j,c,schnitt;
  c=0;
  for(i=0;i<siglength*5;i++)
  {
    if(MaxCorrs[i]>cutoff)
    {
      for(j=i+5;j<siglength*5;j++)
      {
        if(MaxCorrs[j]>cutoff)
        {
          schnitt=Schnitt(Groups[i],Groups[j]);
          if(schnitt>mincov/2 && Groupsize(Groups[i])>schnitt*2 && Groupsize(Groups[j])>schnitt*2)
          {
            Inplace_Schnitt_Group(Groups[c],Groups[i],Groups[j]);
            GrOne(LocalCoverage[c/5]);
            MaxCorrs[c]=cutoff+1.0;
            c++;
            if(c==siglength*5)return;
          }
        }
      }
    }
  }
  for(i=c;i<siglength*5;i++)MaxCorrs[i]=0.0;
  return;
}

// Hier werden die Groups mit den korrelierenden Groups refined und alle Datenstrukturen dafür gefüllt.
void Group_Refinement(double *MaxCorrs, double cutoff,int anfang, int ende, int mincov, int maxclique, double greedy)
{
  int i;

  for(i=0;i<siglength*5;i++)
  {
    Sizes[i]=0;
    Cutoffs[i]=0;
    C_Groups[i]=NULL;
    C_Coverage[i]=NULL;
    Cliques[i]=NULL;
    Drop_Off[i]=1000.0;

    if(MaxCorrs[i]>cutoff)
    {
      Cliques[i]=Cliquer(anfang, ende, mincov, maxclique, greedy, i);
      while(Cliques[i][Sizes[i]]>0){Sizes[i]++;}

      if(Sizes[i]>5)
      {
        #ifdef PRINT
        printf("%d\n",i);
        #endif
        //for(int k=0;k<Sizes[i];k++)printf("%d ",Cliques[i][k] );
        //printf("\n");
        Cutoffs[i]=BestCutoff(30, signumber, Sizes[i], 0.70, 0.05);
        Cutoffs[i]=KorrMaxCutoff(i,Cutoffs[i]);
        Cutoffs[i]=Dropoff_Cutoff(i, 0);
        C_Groups[i]=CliqueGroup(Cliques[i], Cutoffs[i]);
        //C_Core[i]=CoreGroup(Cliques[i], Core_Cutoff(Groupsizearray[i]) );
        C_Coverage[i]=CliqueCoverage(Cliques[i], Cutoffs[i]);
        #ifdef PRINT
        printf("gr->C_gr %d->%d\n",Groupsize(Groups[i]),Groupsize(C_Groups[i]) );
        #endif

        // HIER EIN TEST:
/*        GroupPrecision(Groups[i]);
        GroupPrecision(C_Groups[i]);
        printf("Sizes_Dropoff_Cutoff:%d %d -> %f\n",Sizes[i],Cutoffs[i],Drop_Off[i]);
        Cutoffs[i]=Sizes_Dropoff_Cutoff(i, 1, 5);
        C_Groups[i]=CliqueGroup(Cliques[i], Cutoffs[i]);
        printf("Nachher:%d %d -> %f\n",Sizes[i],Cutoffs[i],Drop_Off[i]);*/

        GroupPrecision(Groups[i]);
        GroupPrecision(C_Groups[i]);
        #ifdef PRINT
        printf("\n");
        #endif

      }
      else
      {
        MaxCorrs[i]=0.0;
      }
    }
  }
}

// Eine parallelisierte Version des GroupRefinement

double *ParaCorrs;
void *Hilfs_Group_Refinement(void *x)
{

  int * args;
  args=((int *) x);

  int anfang=(*(args+0));
  int ende=(*(args+1)); 
  int mincov=(*(args+2));
  int maxclique=(*(args+3)); 
  double cutoff=(*(args+4)); 
  double greedy=(*(args+5)); 
  int thread=(*(args+6)); 
  int NTHREADS=(*(args+7)); 

  int i;

  for(i=0;i<siglength*5;i++)
  {
    if(i%NTHREADS==thread)
    {
      Sizes[i]=0;
      Cutoffs[i]=0;
      C_Groups[i]=NULL;
      C_Coverage[i]=NULL;
      Cliques[i]=NULL;
      Drop_Off[i]=1000.0;

      if(ParaCorrs[i]>cutoff)
      {
        Cliques[i]=Cliquer(anfang, ende, mincov, maxclique, greedy, i);
        while(Cliques[i][Sizes[i]]>0){Sizes[i]++;}

        if(Sizes[i]>5)
        {
          #ifdef PRINT
          printf("%d\n",i);
          #endif
          //for(int k=0;k<Sizes[i];k++)printf("%d ",Cliques[i][k] );
          //printf("\n");
          Cutoffs[i]=BestCutoff(30, signumber, Sizes[i], 0.70, 0.05);
          Cutoffs[i]=KorrMaxCutoff(i,Cutoffs[i]);
          Cutoffs[i]=Dropoff_Cutoff(i, 0);
          C_Groups[i]=CliqueGroup(Cliques[i], Cutoffs[i]);
          //C_Core[i]=CoreGroup(Cliques[i], Core_Cutoff(Groupsizearray[i]) );
          C_Coverage[i]=CliqueCoverage(Cliques[i], Cutoffs[i]);
          #ifdef PRINT
          printf("gr->C_gr %d->%d\n",Groupsize(Groups[i]),Groupsize(C_Groups[i]) );
          #endif

          // HIER EIN TEST:
  /*        GroupPrecision(Groups[i]);
          GroupPrecision(C_Groups[i]);
          printf("Sizes_Dropoff_Cutoff:%d %d -> %f\n",Sizes[i],Cutoffs[i],Drop_Off[i]);
          Cutoffs[i]=Sizes_Dropoff_Cutoff(i, 1, 5);
          C_Groups[i]=CliqueGroup(Cliques[i], Cutoffs[i]);
          printf("Nachher:%d %d -> %f\n",Sizes[i],Cutoffs[i],Drop_Off[i]);*/

          GroupPrecision(Groups[i]);
          GroupPrecision(C_Groups[i]);
          #ifdef PRINT
          printf("\n");
          #endif

        }
        else
        {
          ParaCorrs[i]=0.0;
        }
      }
    }
  }
  return NULL;
}

void Parallel_Group_Refinement(double *MaxCorrs, double cutoff,int anfang, int ende, int mincov, int maxclique, double greedy, int NTHREADS)
{
  time_t start_time=time(NULL);

  pthread_t threads[NTHREADS];
  int *thread_args[NTHREADS];

  int rc,i;

  #ifdef PRINT
  printf("Parallel_Group_Refinement\n");
  #endif
  ParaCorrs=MaxCorrs;

  /* spawn the threads */
  for (i=0; i<NTHREADS; ++i)
  {
    //(int anfang, int ende, int mincov, int maxgroup, double cutoff, int NTHREADS, int thread)
    thread_args[i] = Guarded_Malloc(sizeof(int)*8);
    thread_args[i][0]=anfang;
    thread_args[i][1]=ende;
    thread_args[i][2]=mincov;
    thread_args[i][3]=maxclique;
    thread_args[i][4]=cutoff;
    thread_args[i][5]=greedy;
    thread_args[i][6]=i;
    thread_args[i][7]=NTHREADS;

    #ifdef PRINT
    printf("spawning thread %d\n", i);
    #endif
    rc = pthread_create(&threads[i], NULL, Hilfs_Group_Refinement, (void *) thread_args[i]);
  }

  #ifdef PRINT
  printf("threads done\n");fflush(stdout);
  #endif

  /* wait for threads to finish */
  for (i=0; i<NTHREADS; ++i) {
    rc = pthread_join(threads[i], NULL);
  }
  #ifdef PRINT
  printf("pthread_join\n");fflush(stdout);
  #endif

  time_t time_spent=time(NULL)-start_time;
  #ifdef PRINT
  printf("%ld sec.\n", time_spent);  
  #endif

}

int Unterteilungskomprimierung(int *Unterteilung)
{
  //printf("Unterteilungskomprimierung\n");
  //fflush(stdout);
  int i;
  int max=0;
  for(i=0;i<signumber;i++){if(max<Unterteilung[i])max=Unterteilung[i];}
  int *Replace=Guarded_Malloc(sizeof(int)*(max+1)); 
  for(i=0;i<max+1;i++)Replace[i]=-1; 
  max=0;
  for(i=0;i<signumber;i++)
  {
    if(Unterteilung[i]>-1)
    {
      if(Replace[Unterteilung[i]]<0){Replace[Unterteilung[i]]=max;max++;}
      Unterteilung[i]=Replace[Unterteilung[i]];
    }
  }
  free(Replace);
  return max;
}

int* UnterteilungsKomplettierung(int *Unterteilung)
{
  //printf("Unterteilungskomprimierung\n");
  //fflush(stdout);
  int *NewUnterteilung=Guarded_Malloc(realsigno*sizeof(int));
  int i,j;
  j=0;
  for(i=0;i<realsigno;i++)
  {
    if(Ausgelassen[i]==1)
    {
      NewUnterteilung[i]=Unterteilung[j];
      j++;
    }
    if(Ausgelassen[i]==-1)
    {
      NewUnterteilung[i]=-1;
    }
  }
  return NewUnterteilung;
}

void Matrixprint(double **Matrix, int x, int y)
{
  int i,j;
  for(i=0;i<x;i++)
  {
    for(j=0;j<y;j++)
    {
      printf("%.2f ",Matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

// Hier berechne ich Wahrscheinlichkeiten für die Verbindungen statt einfach welche auszugeben
void Verbindungen(int *Unterteilung1, int *Unterteilung2, int *Unterteilung3)
{
  int i,j,k;
  int max1=Unterteilungskomprimierung(Unterteilung1);
  int max2=Unterteilungskomprimierung(Unterteilung2);
  int max3=Unterteilungskomprimierung(Unterteilung3);

  #ifdef PRINT
  printf("Maxes: %d %d %d\n",max1,max2,max3);
  for(i=0;i<signumber;i++){printf("%d ",Unterteilung1[i]);}
  printf("\n\n");
  for(i=0;i<signumber;i++){printf("%d ",Unterteilung2[i]);}
  printf("\n\n");
  for(i=0;i<signumber;i++){printf("%d ",Unterteilung3[i]);}
  printf("\n\n");
  #endif

  double **Matrix12=Guarded_Malloc(sizeof(double*)*max1);
  for(i=0;i<max1;i++)
  {
    Matrix12[i]=Guarded_Malloc(sizeof(double)*max2);
    for(j=0;j<max2;j++)Matrix12[i][j]=0.0;
  }

  double **Matrix23=Guarded_Malloc(sizeof(double*)*max2);
  for(i=0;i<max2;i++)
  {
    Matrix23[i]=Guarded_Malloc(sizeof(double)*max3);    
    for(j=0;j<max3;j++)Matrix23[i][j]=0.0;
  }

  double **Matrix13=Guarded_Malloc(sizeof(double*)*max1);
  for(i=0;i<max1;i++)
  {
    Matrix13[i]=Guarded_Malloc(sizeof(double)*max3);    
    for(j=0;j<max3;j++)Matrix13[i][j]=0.0;
  }

  double *Card1=Guarded_Malloc(sizeof(double)*max1);
  for(i=0;i<max1;i++)Card1[i]=0.0;

  double *Card2in=Guarded_Malloc(sizeof(double)*max2);
  for(i=0;i<max2;i++)Card2in[i]=0.0;
  double *Card2out=Guarded_Malloc(sizeof(double)*max2);
  for(i=0;i<max2;i++)Card2out[i]=0.0;    

  double *Card3=Guarded_Malloc(sizeof(double)*max3);
  for(i=0;i<max3;i++)Card3[i]=0.0;

  for(i=0;i<signumber;i++)
  {
    if(Unterteilung1[i]>-1)
    {
      Card1[Unterteilung1[i]]+=1.0;
      if(Unterteilung2[i]>-1)
      {
        Card2in[Unterteilung2[i]]+=1.0;
        Matrix12[Unterteilung1[i]][Unterteilung2[i]]+=1.0;
      }
    }
    if(Unterteilung3[i]>-1)
    {
      Card3[Unterteilung3[i]]+=1.0;
      if(Unterteilung2[i]>-1)
      {
        Card2out[Unterteilung2[i]]+=1.0;
        Matrix23[Unterteilung2[i]][Unterteilung3[i]]+=1.0;
      }
    }

  }
/*  printf("Die Matrizen\n");fflush(stdout);
  Matrixprint(Matrix12, max1, max2);
  Matrixprint(Matrix23, max2, max3);
  Matrixprint(Matrix13, max1, max3);*/

  // Normalisierte Vorwärtsverbindungen sind Matrix12[i][k]/Card1[i] und Matrix23[k][j]/Card2out[k]
  double sum;
  double back;
  for(i=0;i<max1;i++) //Alle linken Gruppen
  {
    for(j=0;j<max3;j++) //Alle rechten Gruppen
    {
      sum=0.0;
      back=0.0;
      for(k=0;k<max2;k++) //Alle möglichen Verbindungen zwischen ihnen
      {
        if(Card1[i]>0.5 && Card2out[k]>0.5)
        {
          sum+=(Matrix12[i][k]/Card1[i])*(Matrix23[k][j]/Card2out[k]);
        }
        if(Card2in[k]>0.5 && Card3[j]>0.5) //Der Rückweg? 
        {
          back+=(Matrix12[i][k]/Card2in[k])*(Matrix23[k][j]/Card3[j]);
        }
      }
      Matrix13[i][j]=sum*back; // Oder max? Oder 1.0-(1.0-sum)(1.0-back)?
    }
  }
  //Ausgabe:
  double max;
  int maxj;
  int count=0;
  int count2=0;
  int count3=0;
  int count4=0;
  int count5=0;
  int falsep=0;
  int falsep2=0;
  int falsep3=0;
  int falsep4=0;
  int falsep5=0;  
  for(i=0;i<max1;i++) //Alle linken Gruppen
  {
    maxj=0;
    max=0.0;
    for(j=0;j<max3;j++) //Alle rechten Gruppen
    {
      if(Matrix13[i][j]>max)
      {
        max=Matrix13[i][j];
        maxj=j;
      }
    }
    #ifdef PRINT
    printf("%d -> %d : %f\n",i,maxj,max );
    #endif
    if(i==maxj)count++;
    if(max>0.25 && i==maxj)count2++;
    if(max>0.50 && i==maxj)count3++;
    if(max>0.75 && i==maxj)count4++;
    if(max>0.90 && i==maxj)count5++;

    if(i!=maxj)falsep++;
    if(max>0.25 && i!=maxj)falsep2++;
    if(max>0.50 && i!=maxj)falsep3++;
    if(max>0.75 && i!=maxj)falsep4++;
    if(max>0.90 && i!=maxj)falsep5++;    
  }  
  #ifdef PRINT
  printf("%d, %d, %d, %d, %d korrekt mit > 0.0, 0.25, 0.50, 0.75, 0.90. \n",count,count2,count3,count4,count5 );
  printf("%d, %d, %d, %d, %d false positives > 0.0, 0.25, 0.50, 0.75, 0.90. \n",falsep,falsep2,falsep3,falsep4,falsep5 );
  #endif
}


// Hier berechne ich Wahrscheinlichkeiten für die Verbindungen statt einfach welche auszugeben
double** Multi_Verbindungen(int **Unterteilungen, int unt_no)
{
  int i,j,k;
  int max[unt_no];
  for(i=0;i<unt_no;i++){max[i]=Unterteilungskomprimierung(Unterteilungen[i]);}
  int maxmax=0;
  for(i=0;i<unt_no;i++){if(maxmax<max[i])maxmax=max[i];}

  double **Matrices[unt_no];
  for(j=0;j<unt_no-1;j++)
  {
    Matrices[j]=Guarded_Malloc(sizeof(double*)*max[j]);
    for(i=0;i<max[j];i++)
    {
      Matrices[j][i]=Guarded_Malloc(sizeof(double)*max[j+1]);
      for(k=0;k<max[j+1];k++)Matrices[j][i][k]=0.0;
    }    
  }

  double **Matrix[2];
  Matrix[0]=Guarded_Malloc(sizeof(double*)*maxmax); //Für die Zwischenergebnisse
  Matrix[1]=Guarded_Malloc(sizeof(double*)*maxmax); //Für die Zwischenergebnisse
  for(i=0;i<maxmax;i++)
  {
    Matrix[0][i]=Guarded_Malloc(sizeof(double)*maxmax);
    Matrix[1][i]=Guarded_Malloc(sizeof(double)*maxmax);
    for(k=0;k<maxmax;k++)
    {
      Matrix[1][i][k]=0.0;  
      Matrix[0][i][k]=0.0; 
    }   
  }

  double *Card_out[unt_no-1];
  for(j=0;j<unt_no-1;j++)
  {
    Card_out[j]=Guarded_Malloc(sizeof(double)*max[j]);
    for(i=0;i<max[j];i++)
    {
      Card_out[j][i]=0.0;
    }
  }

  //for(j=0;j<unt_no;j++)printf("%d\n",max[j] );
  #ifdef PRINT
  printf("Das Allozieren\n");
  fflush(stdout);
  #endif 

  for(i=0;i<signumber;i++)
  {
    //printf("\n%d - ",i );fflush(stdout);
    for(j=0;j<unt_no-1;j++)
    {
      //printf("%d ",Unterteilungen[j][i] );fflush(stdout);
      //if( (Unterteilungen[j][i]>-1 || j==0) && (Unterteilungen[j+1][i]>-1 || j+1==unt_no-1) ) //Keine Ahnung was das sollte
      if(Unterteilungen[j][i]>-1 && Unterteilungen[j+1][i]>-1)
      {
        Card_out[j][Unterteilungen[j][i]]+=1.0;
        Matrices[j][Unterteilungen[j][i]][Unterteilungen[j+1][i]]+=1.0;
      }
    }
  }

  #ifdef PRINT
  printf("Das Füllen\n");
  fflush(stdout);
  #endif


  //Jetzt die Berechnung: 
  double sum;
  int l=0;  //Erst wird die erste Matrix relativiert in die Zwischenresultmatrix geschrieben
  for(i=0;i<max[l];i++) 
  {
    for(j=0;j<max[l+1];j++) 
    {  
      if(Card_out[0][i]>0.5)Matrix[1][i][j]=Matrices[0][i][j]/Card_out[0][i];
    }
  }
  for(l=0;l<unt_no-2;l++)
  {
    for(i=0;i<max[l];i++) //Alle linken Gruppen
    {
      for(j=0;j<max[l+2];j++) //Alle rechten Gruppen
      {
        sum=0.0;
        for(k=0;k<max[l+1];k++) //Alle möglichen Verbindungen zwischen ihnen
        {
          if(Card_out[l][i]>0.5 && Card_out[l+1][k]>0.5)
          {
            sum+=Matrix[(l+1)%2][i][k]*(Matrices[l+1][k][j]/Card_out[l+1][k]);
          }
        }
        Matrix[l%2][i][j]=sum; // Oder max? Oder 1.0-(1.0-sum)(1.0-back)?
      }
    }
  }

  #ifdef PRINT
  printf("Die Berechnung\n");
  fflush(stdout);
  #endif

  //Free everything
  for(j=0;j<unt_no-1;j++)
  {
    free(Card_out[j]);
  }
  for(j=0;j<unt_no-1;j++)
  {
    for(i=0;i<max[j];i++)
    {
      free(Matrices[j][i]);
    }    
    free(Matrices[j]);
  }
  free(Matrix[(unt_no)%2]);

  #ifdef PRINT
  printf("Das Freeen\n");
  fflush(stdout);
  #endif

  //Das Ergebnis:
/*  for(i=0;i<max[0];i++)
  {
    for(j=0;j<max[unt_no-1];j++)
    {
      printf("%f ",Matrix[(unt_no-1)%2][i][j]);
    }
    printf("\n");
  }*/

  //Except the result:
  return Matrix[(unt_no-1)%2];
}




void Back_And_Forth(int **Unterteilungen, int unt_no)
{
  int i,j;
  double **Matrix_Hin=Multi_Verbindungen(Unterteilungen, unt_no);
  int **Unterteilungen2=Guarded_Malloc(sizeof(int*)*unt_no);
  for(i=0;i<unt_no;i++)Unterteilungen2[unt_no-i-1]=Unterteilungen[i];
  double **Matrix_Back=Multi_Verbindungen(Unterteilungen2, unt_no);
  int max[unt_no];
  for(i=0;i<unt_no;i++){max[i]=Unterteilungskomprimierung(Unterteilungen[i]);}

  //Was ist die Verbindung von i<max[0] zu j<max[unt_no-1]?: Matrix_Hin[i][j]*Matrix_Back[j][i];  


  //Ausgabe:
  double maxx;
  int maxj;
  int count=0;
  int count2=0;
  int count3=0;
  int count4=0;
  int count5=0;
  int falsep=0;
  int falsep2=0;
  int falsep3=0;
  int falsep4=0;
  int falsep5=0;  
  for(i=0;i<max[0];i++) //Alle linken Gruppen
  {
    maxj=0;
    maxx=0.0;
    for(j=0;j<max[unt_no-1];j++) //Alle rechten Gruppen
    {
      if(Matrix_Hin[i][j]*Matrix_Back[j][i]>maxx)
      {
        maxx=Matrix_Hin[i][j]*Matrix_Back[j][i];
        maxj=j;
      }
    }

    printf("%d -> %d : %f\n",i,maxj,maxx );

    if(i==maxj)count++;
    if(maxx>0.25 && i==maxj)count2++;
    if(maxx>0.50 && i==maxj)count3++;
    if(maxx>0.75 && i==maxj)count4++;
    if(maxx>0.90 && i==maxj)count5++;

    if(i!=maxj)falsep++;
    if(maxx>0.25 && i!=maxj)falsep2++;
    if(maxx>0.50 && i!=maxj)falsep3++;
    if(maxx>0.75 && i!=maxj)falsep4++;
    if(maxx>0.90 && i!=maxj)falsep5++;    
  }  
  printf("%d, %d, %d, %d, %d korrekt mit > 0.0, 0.25, 0.50, 0.75, 0.90. \n",count,count2,count3,count4,count5 );
  printf("%d, %d, %d, %d, %d false positives > 0.0, 0.25, 0.50, 0.75, 0.90. \n",falsep,falsep2,falsep3,falsep4,falsep5 );

}



// Hier unterteile ich die Sigs nach und nach nach den besten C_Groups:
void Subdivision(double *MaxCorrs, double cutoff,int sizecutoff,int anfang, int ende, int mingroup)
{
  //Sortieren der C_Groups nach Drop_Off, Corr average und Size:
  int anzahl=0;
  int *I=Guarded_Malloc(sizeof(int)*siglength*5); //Indizes
  int i,j,k;
  for(i=0;i<siglength*5;i++)
  {
    if(MaxCorrs[i]>cutoff && Sizes[i]>sizecutoff)
    {
      I[anzahl]=i;
      anzahl++;
    }
  }
  #ifdef PRINT
  printf("Es gibt %d Groups über den Cutoffs.\n",anzahl );
  #endif
  for(i=0;i<anzahl;i++)
  {
    for(j=i+1;j<anzahl;j++)
    {
      if(Drop_Off[I[i]]>Drop_Off[I[j]]){k=I[i];I[i]=I[j];I[j]=k;}
      else if(Drop_Off[I[i]]==Drop_Off[I[j]])
      {
        if(Sizes[I[i]]<Sizes[I[j]]){k=I[i];I[i]=I[j];I[j]=k;}
        else if(Sizes[I[i]]==Sizes[I[j]])
        {
          if(MaxCorrs[I[i]]<MaxCorrs[I[j]]){k=I[i];I[i]=I[j];I[j]=k;}
        }
      }
    }
  }
  #ifdef PRINT
  printf("Sortiert.\n");
  #endif

  //Sukzessives Unterteilen:
  int *Unterteilung=Guarded_Malloc(sizeof(int)*signumber); //Unterteilung
  for(i=0;i<signumber;i++)Unterteilung[i]=0;
  int number=1;
  int number2=1;
  int drinne,draus;
  for(i=0;i<anzahl;i++) //Durch die C_Groups
  {
    for(k=0;k<number;k++) // Alle bisherigen Aufteilungen werden auf Aufteilbarkeit untersucht:
    {
      drinne=0;
      draus=0;
      for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
      {
        if(Unterteilung[j]==k) //Ist in Aufteilung k
        {
          //printf("%d %d %d\n",i,I[i],j);
          //fflush(stdout);
          if(GrElement(C_Groups[I[i]],j)) {drinne+=1;}
          else {draus+=1;}
        }
      }
      //printf("gr %d -> drinne %d draus %d\n",Groupsize(C_Groups[I[i]]), drinne,draus);
      if(drinne>mingroup && draus>mingroup) //Wenn die Auteilung nicht mit bisherigen clasht wird aufgeteilt:
      {
        for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
        {
          if(Unterteilung[j]==k) //Ist in Aufteilung k
          {
            //printf("- %d %d %d\n",i,I[i],j);
            //fflush(stdout);
            if(GrElement(C_Groups[I[i]],j)) {Unterteilung[j]=number2;}
            else {Unterteilung[j]=number2+1;}
          }
        }
        number2+=2; 
      }
    }
    number=number2;
    #ifdef PRINT
    printf("Nach %d von %d: number der Unterteilungen %d.\n",i,anzahl,number );
    fflush(stdout);
    #endif
  }  

  #ifdef PRINT
  printf("Unterteilt.\n");
  fflush(stdout);
  #endif

  return; 

  //Assessment: 
  for(i=0;i<signumber;i++)printf("%d -> %d \n",i,Unterteilung[i]);

  int *Counter=Guarded_Malloc(sizeof(int)*(number+1));
  int *Comparer=Guarded_Malloc(sizeof(int)*(30));
  for(i=0;i<number+1;i++)Counter[i]=0;
  for(i=0;i<signumber;i++)Counter[Unterteilung[i]]++;

  printf("Unterteilungsgrößen:\n");
  for(i=0;i<number+1;i++)printf("%d -> %d\n",i,Counter[i]);
  printf("\n");

  for(i=0;i<signumber/30;i++)
  {
    for(j=0;j<30;j++)Comparer[j]=0;
    for(j=0;j<30;j++)
    {
      for(k=0;k<30;k++)
      {
        if(Unterteilung[i*30+j]==Unterteilung[i*30+k])Comparer[j]++;
      }
    }
    draus=-1;
    drinne=0;
    for(j=0;j<30;j++){if(drinne<Comparer[j])drinne=Comparer[j];}
    for(j=0;j<30;j++)if(Comparer[j]==drinne){draus=Unterteilung[i*30+j];k=j;} //Das ist die Hauptbeschreibung
    if(draus==-1)printf("Gruppe %d: Es gibt keine Beschreibung.\n",i);
    else {printf("Gruppe %d: Drin %d / %d\n",i,Comparer[k], Counter[Unterteilung[i*30+k]]);}
  }

  // Hier verbinde ich halbe Gruppen nach größter Übereinstimmung:
  int **Matrix=Guarded_Malloc(sizeof(int*)*100);
  for(i=0;i<100;i++)Matrix[i]=Guarded_Malloc(sizeof(int)*100);

  int maxj,maxi,maxshare,sharecount;
  int ii,jj;
  int correct_connect=0;
  for(i=0;i<signumber/30;i++)
  {
    maxj=-1;
    maxshare=0;
    for(j=0;j<signumber/30;j++)
    {
      sharecount=0;
      for(ii=0;ii<15;ii++)
      {
        for(jj=0;jj<15;jj++)
        {
          if(Unterteilung[i*30+ii]==Unterteilung[j*30+15+jj])sharecount++;
        }
      }
      if(sharecount>maxshare)
      {
        maxshare=sharecount;
        maxj=j;
      }
      Matrix[i][j]=sharecount;
    }
    printf("%d -> %d\n",i,maxj);
    if(i==maxj)correct_connect++;
  }
  printf("%d korrekte Verbindungen.\n",correct_connect);

  //Bessere Auflösung:
  int *Resolution=Guarded_Malloc(sizeof(int)*100);
  for(i=0;i<100;i++)Resolution[i]=-1;

  for(ii=0;ii<100;ii++) //Der i-th beste Eintrag:
  {
    maxi=-1;
    maxj=-1;
    maxshare=0;
    for(j=0;j<100;j++)
    {
      for(i=99;i>-1;i--)
      {
        if(Matrix[j][i]>maxshare)
        {
          maxshare=Matrix[j][i];
          maxi=i;
          maxj=j;
        }
      }
    }
    if(maxi>-1 && maxj>-1)
    {
      Resolution[maxj]=maxi;
      for(i=0;i<100;i++)
      {
        Matrix[maxj][i]=0;
        Matrix[i][maxi]=0;
      }
    }
  }

  for(i=0;i<100;i++)printf("%d -> %d\n",i,Resolution[i] );
  correct_connect=0;
  for(i=0;i<100;i++){if(i==Resolution[i])correct_connect++;}
  printf("%d korrekte Verbindungen.\n",correct_connect);

  return;
}

// Hier kommt die Implementierung des KMeans-Zeugs. 
// Erstmal herausfinden, welche Vars man für eine Aufteilung verwenden soll:
int *Relative_Vars(int *Unterteilung, int u_no, double *MaxCorrs, double cutoff, int mingroup)
{
  int i,j,count;
  unsigned long *U_Group=GrInitialize();
  int *SelectedVars=Guarded_Malloc(sizeof(int)*siglength*5);
  count=0;
  for(i=0;i<siglength*5;i++)
  {
    if(MaxCorrs[i]>cutoff){SelectedVars[i]=1;count++;}
    else{SelectedVars[i]=0;}
  }
  #ifdef PRINT
  printf("Vars %d\n",count);
  #endif
  for(i=0;i<signumber;i++)if(Unterteilung[i]==u_no)GrAdd(U_Group,i);

  #ifdef PRINT
  printf("Size of Subgroup: %d\n",Groupsize(U_Group));
  #endif

  for(i=0;i<siglength*5;i++) //Die nicht genug in der Unterteilung sind
  {
    if(SelectedVars[i])
    {
      if(Schnitt(U_Group,Groups[i])<mingroup){SelectedVars[i]=0;count--;}
    }
  }
  #ifdef PRINT
  printf("Vars %d\n",count);
  #endif

  //Jetzt werden die ausgesiebt, die nicht innerhalb von u_no korrelieren:
  double Z;
  count=0;
  for(i=0;i<siglength*5;i++)
  {
    if(SelectedVars[i])
    {
      for(j=i+100;j<siglength*5;j++)
      {
        if(SelectedVars[j])
        {
          Z=Relative_Group_Significance(Groups[j],Groups[i],U_Group);
          //printf("%f %d %d\n",Z,i,j );
          if(Z>cutoff)
          {
            SelectedVars[i]=2;
            SelectedVars[j]=2;
          }
        }
      }
    }
  }

  free(U_Group);
  for(i=0;i<siglength*5;i++){if(SelectedVars[i]==2)count++;}
  #ifdef PRINT
  printf("Vars %d\n",count);
  #endif
  int *Vars=Guarded_Malloc(sizeof(int)*(count+1));
  Vars[count]=-1;
  count=0;
  for(i=0;i<siglength*5;i++){if(SelectedVars[i]==2){Vars[count]=i;count++;}}
  free(SelectedVars);
  #ifdef PRINT
  for(j=0;j<count;j++)printf("%d ",Vars[j] );
  printf("\n");
  #endif
  return Vars;
}


// Dann ein Clustering, das mehr oder weniger Kmeans ist:
void Kmeans1(int *Unterteilung, int u_no, int *Vars, int mingroup)
{
  int anzahl=0;
  int varzahl=0;
  while(Vars[varzahl]!=-1)varzahl++;
  #ifdef PRINT
  printf("Varanzahl %d\n",varzahl);
  #endif
  int old_sc=sc;
  int *I=Guarded_Malloc(sizeof(int)*signumber);
  int i,j;

  for(i=0;i<signumber;i++) // anzahl und Indices der sigs
  {
    if(Unterteilung[i]==u_no)
    {
      I[anzahl]=i;
      anzahl++;
    }
  }

  sc=(varzahl/64)+1;
  #ifdef PRINT
  printf("New sc %d, varzahl %d, anzahl= %d / %d\n",sc,varzahl,anzahl,signumber);
  for(j=0;j<varzahl;j++)printf("%d ",Vars[j] );
  printf("\n");
  #endif

  unsigned long **VarSigs=Guarded_Malloc(sizeof(unsigned long*)*anzahl);
  for(i=0;i<anzahl;i++) //VarSigs als Bitvektoren erstellen
  {
    VarSigs[i]=GrInitialize();
    for(j=0;j<varzahl;j++) //Durch die Vars
    {
      if(GrElement(Groups[Vars[j]],I[i])){GrAdd(VarSigs[i],j);}
    }
  }


  //Dann das Clustern:
  int *Clusternumber=Guarded_Malloc(sizeof(int)*anzahl);
  int *Clustersize=Guarded_Malloc(sizeof(int)*anzahl);
  for(i=0;i<anzahl;i++)Clustersize[i]=0;
  int best_j,best_score,score;
  for(i=0;i<anzahl;i++) //Verteilen auf die anderen Sigs
  {
    best_score=0;
    best_j=0;
    for(j=0;j<anzahl;j++)
    {
      score=GrMatch(VarSigs[j], VarSigs[i]);

      if(score>best_score && i!=j)
      {
        //printf("%d %d %d\n",i,j,score );
        best_score=score;
        best_j=j;
      }
    }
    Clusternumber[i]=best_j;
    Clustersize[best_j]++;
  }

  int min;
  for(min=1;min<mingroup;min++)
  {
    for(i=0;i<anzahl;i++) //Verteilen der Kleinen auf die Großen:
    {
      if(Clustersize[Clusternumber[i]]<min) //i in einem zu kleinen Cluster
      {
        best_score=0;
        best_j=0;
        for(j=0;j<anzahl;j++) // Durch die Cluster
        {
          if(Clustersize[j]>=min)  // Cluster groß genug
          {
            score=GrMatch(VarSigs[j], VarSigs[i]);
            if(score>best_score && i!=j)
            {
              best_score=score;
              best_j=j;
            }
          }
        }
        Clusternumber[i]=best_j;
        Clustersize[best_j]++;
      }
    }
  }


  free(I);
  free(Clustersize);
  for(i=0;i<anzahl;i++)free(VarSigs[i]);
  free(VarSigs);

  sc=old_sc;

  //Ergebnis in die Unterteilung eintragen:
  int max_u=0;
  for(i=0;i<signumber;i++){if(Unterteilung[i]>max_u)max_u=Unterteilung[i];}
  for(i=0;i<anzahl;i++)Unterteilung[I[i]]=Clusternumber[i]+max_u+1;
  free(Clusternumber);

}

// Dann ein Clustering, das mehr oder weniger Kmeans ist: Erfolg der Aufteilung wird returned:
int Kmeans(int *Unterteilung, int u_no, int *Vars, int mingroup)
{
  int anzahl=0;
  int varzahl=0;
  while(Vars[varzahl]!=-1)varzahl++;
  #ifdef PRINT
  printf("Varanzahl %d\n",varzahl);
  #endif
  int old_sc=sc;
  int *I=Guarded_Malloc(sizeof(int)*signumber);
  int i,j,k;

  for(i=0;i<signumber;i++) // anzahl und Indices der sigs
  {
    if(Unterteilung[i]==u_no)
    {
      I[anzahl]=i;
      anzahl++;
    }
  }

  sc=(varzahl/64)+1;
  #ifdef PRINT
  printf("New sc %d, varzahl %d, anzahl= %d / %d\n",sc,varzahl,anzahl,signumber);
  for(j=0;j<varzahl;j++)printf("%d ",Vars[j] );
  printf("\n");
  #endif

  unsigned long **VarSigs=Guarded_Malloc(sizeof(unsigned long*)*anzahl);
  for(i=0;i<anzahl;i++) //VarSigs als Bitvektoren erstellen
  {
    VarSigs[i]=GrInitialize();
    for(j=0;j<varzahl;j++) //Durch die Vars
    {
      if(GrElement(Groups[Vars[j]],I[i])){GrAdd(VarSigs[i],j);}
    }
  }

  unsigned long **Centroids=Guarded_Malloc(sizeof(unsigned long*)*anzahl);
  for(i=0;i<anzahl;i++) //VarSigs als Bitvektoren erstellen
  {
    Centroids[i]=GrInitialize();
  }

  //Dann das Centroiden:
  int *Clusternumber=Guarded_Malloc(sizeof(int)*anzahl);
  int *Clustersize=Guarded_Malloc(sizeof(int)*anzahl);
  for(i=0;i<anzahl;i++)Clustersize[i]=0;
  int best_j[5];
  int best_score[5];
  int score;
  int l,s;
  for(i=0;i<anzahl;i++) //Verteilen auf die anderen Sigs
  {
    for(k=0;k<5;k++)
    {
      best_score[k]=0;
      best_j[k]=0;
    }

    for(j=0;j<anzahl;j++)
    {
      score=GrMatch(VarSigs[j], VarSigs[i]);
      for(k=0;k<5;k++)
      {
        for(l=k+1;l<5;l++)
        {
          if(best_score[l]<best_score[k])
          {
            s=best_score[l];
            best_score[l]=best_score[k];
            best_score[k]=s;
            s=best_j[l];
            best_j[l]=best_j[k];
            best_j[k]=s;
          }
        }
      }
      if(score>best_score[0])
      {
        //printf("%d %d %d\n",i,j,score );
        best_score[0]=score;
        best_j[0]=j;
      }
    }
    //printf("%d = ",i );
    //for(k=0;k<5;k++)printf("%d-%d ",best_score[k],best_j[k]);
    //printf("\n");    

    //Hier die Erstellung des Centroids:
    for(j=0;j<varzahl;j++) //Durch die Vars
    {
      s=0;
      for(k=0;k<5;k++)
      {
        if(GrElement(Groups[Vars[j]],I[best_j[k]])){s++;}
      }
      if(s>2){GrAdd(Centroids[i],j);}
    }
  }

  //Dann das Clustern:
  for(i=0;i<anzahl;i++) //Verteilen auf die anderen Sigs
  {
    best_score[0]=0;
    best_j[0]=0;
    for(j=0;j<anzahl;j++)
    {
      score=GrMatch(Centroids[j], VarSigs[i]);

      if(score>best_score[0] && i!=j)
      {
        //printf("%d %d %d\n",i,j,score );
        best_score[0]=score;
        best_j[0]=j;
      }
    }
    Clusternumber[i]=best_j[0];
    Clustersize[best_j[0]]++;
  }


  int min;
  for(min=2;min<mingroup;min++)
  {
    for(i=0;i<anzahl;i++) //Verteilen der Kleinen auf die Großen:
    {
      if(Clustersize[Clusternumber[i]]<=min) //i in einem zu kleinen Cluster
      {
        best_score[0]=0;
        best_j[0]=0;
        for(j=0;j<anzahl;j++) // Durch die Cluster
        {
          if(Clustersize[j]>=min && Clusternumber[i]!=j)  // Cluster groß genug
          {
            score=GrMatch(Centroids[j], VarSigs[i]);
            if(score>best_score[0] && i!=j)
            {
              best_score[0]=score;
              best_j[0]=j;
            }
          }
        }
        Clustersize[Clusternumber[i]]--;
        Clusternumber[i]=best_j[0];
        Clustersize[best_j[0]]++;
      }
    }
    //printf("\n%d\n",min);
    //for(j=0;j<10;j++)printf("%d %d %d\n",j,Clusternumber[j],Clustersize[Clusternumber[j]] );
    //fflush(stdout);
  }


  // Mal die Ergebnisse ausgeben:
/*  varzahl=40;
  for(j=0;j<anzahl;j++)
  {
    if(Clustersize[j]>mingroup)
    {
      printf("\n\nCluster %d \n", j);
      for(l=0;l<varzahl;l++)printf("%d",GrElement(Centroids[j],l));
      printf("\n\n");
      for(i=0;i<anzahl;i++)
      {
        if(Clusternumber[i]==j)
        {
          printf("Seq %d, C_no %d, score %d\n",i,Clusternumber[i],GrMatch(Centroids[Clusternumber[i]], VarSigs[i]) );
          for(l=0;l<varzahl;l++)printf("%d",GrElement(VarSigs[j],l));
          printf("\n");
        }
      }
    }
  }
  fflush(stdout);
  exit(0);*/

  int aufgeteilt=0;  //Wurde aufgeteilt oder nicht?
  #ifdef PRINT
  printf("Clustersize: ");
  #endif
  for(i=0;i<anzahl;i++)
  {
    if(Clustersize[i]>0)aufgeteilt++;
    #ifdef PRINT
    printf("%d ",Clustersize[i] );
    #endif
  }
  #ifdef PRINT
  printf("\n");
  #endif


  //TEST:
  #ifdef PRINT
  for(i=0;i<anzahl;i++)
  {
    if(!(I[i]%30==I[i+1]%30))printf("\n");
    printf("%d -> %d(%d)\n",I[i],Clusternumber[i],Clustersize[Clusternumber[i]]);
  }
  #endif

  free(Clustersize);
  for(i=0;i<anzahl;i++)free(VarSigs[i]);
  free(VarSigs);

  sc=old_sc;

  //Ergebnis in die Unterteilung eintragen:
  int max_u=0;
  for(i=0;i<signumber;i++){if(Unterteilung[i]>max_u)max_u=Unterteilung[i];}
  for(i=0;i<anzahl;i++)Unterteilung[I[i]]=Clusternumber[i]+max_u+1;
  free(Clusternumber);
  free(I);

  return aufgeteilt;

}

// Ein Assessment der Unterteilung: 
void Unterteilung_Assessment(int *Unterteilung)
{
  int *U_Count=Guarded_Malloc(sizeof(int)*signumber);
  int *Matrix=Guarded_Malloc(sizeof(int)*30);
  int i,j,k,count,maxi;
  double plus=0;
  double all=0;
  for(i=0;i<signumber;i++)U_Count[i]=0;
  for(i=0;i<signumber;i++)U_Count[Unterteilung[i]]++;
  for(i=0;i<signumber/30;i++)
  {
    for(j=0;j<30;j++)
    {
      Matrix[j]=0;
    }
    for(j=0;j<30;j++)
    {
      for(k=0;k<30;k++)
      {
        if(Unterteilung[i*30+k]==Unterteilung[i*30+j])Matrix[j]++;
      }
    }
    maxi=-1;
    k=0;
    for(j=0;j<30;j++){if(Matrix[j]>maxi){maxi=Matrix[j];k=j;}}
    count=0;
    for(j=0;j<30;j++){if(Unterteilung[i*30+j]==Unterteilung[i*30+k])count++;}  
    //printf("%d/%d\n",count,U_Count[Unterteilung[i*30+k]]);
    plus+=(double)count/(double)U_Count[Unterteilung[i*30+k]];
    all+=1.0;
  }
  printf("U_Acc:%f\n",plus/all );
}

// Der Cutoff, der den Dropoff der C_Groupsizes minimiert:
int Relative_Dropoff_Cutoff(int c_i, int c, int *Unterteilung, int u_no)
{

  unsigned long *Cutoff_Groups[100];
  int i,jj,ii;

  // Die Groups für verschiedene Cutoffs:
  for(i=0;i<Sizes[c_i];i++)
  {
    Cutoff_Groups[i]=GrInitialize();
  }
  for(i=0;i<signumber;i++)
  {
    ii=0; //Der score
    for(jj=0;jj<Sizes[c_i];jj++) //Durch die vars
    {
      if(GrElement(Groups[Cliques[c_i][jj]],i) && Unterteilung[i]==u_no)
      {
        GrAdd(Cutoff_Groups[ii],i);  //In mehr als ii Groups
        ii++;
      }
    }
  }  

  // Minimalen Dropoff:
  double sizes[100];
  for(i=0;i<Sizes[c_i];i++)sizes[i]=(double)Groupsize(Cutoff_Groups[i]);
  int drop_c=intmax(1,c);
  c=drop_c;
  double min_drop=1000000.0; 
  double drop;

  //printf("Sizes:\n");
  //for(i=0;i<Sizes[c_i]-1;i++)printf("%f ",sizes[i]);
  //printf("\n");

  //printf("Drops %d: ",c);
  for(i=c;i<Sizes[c_i]-1;i++)
  {
    if(doublemin((double)signumber-sizes[i],sizes[i])>0)
    {
      drop=((sizes[i-1]-sizes[i+1])/doublemin((double)signumber-sizes[i],sizes[i]));
      //printf("%f ",drop);
      if(drop<min_drop)
      {
        min_drop=drop;
        drop_c=i;
      }
    }
  }
  //printf("\n");

  // Den Wert eintragen.
  Drop_Off[c_i]=min_drop;

  for(i=0;i<Sizes[c_i];i++)
  {
    free(Cutoff_Groups[i]);
  }

  return drop_c;
}

// Hier unterteile ich die Sigs mit relativen Dropoff-Cutoff-C_Groups und benutze nur noch Drop_Off==0
int *Perfect_Subdivision(double *MaxCorrs, double cutoff,int sizecutoff,int anfang, int ende, int mingroup)
{
  //Sortieren der C_Groups nach Drop_Off, Corr average und Size:
  int anzahl=0;
  int *I=Guarded_Malloc(sizeof(int)*siglength*5); //Indizes
  int i,j,k;
  for(i=0;i<siglength*5;i++)
  {
    if(MaxCorrs[i]>cutoff && Sizes[i]>sizecutoff)
    {
      I[anzahl]=i;
      anzahl++;
    }
  }
  #ifdef PRINT
  printf("Es gibt %d Groups über den Cutoffs.\n",anzahl );
  #endif
  for(i=0;i<anzahl;i++)
  {
    for(j=i+1;j<anzahl;j++)
    {
      if(Drop_Off[I[i]]>Drop_Off[I[j]]){k=I[i];I[i]=I[j];I[j]=k;}
      else if(Drop_Off[I[i]]==Drop_Off[I[j]])
      {
        if(Sizes[I[i]]<Sizes[I[j]]){k=I[i];I[i]=I[j];I[j]=k;}
        else if(Sizes[I[i]]==Sizes[I[j]])
        {
          if(MaxCorrs[I[i]]<MaxCorrs[I[j]]){k=I[i];I[i]=I[j];I[j]=k;}
        }
      }
    }
  }
  #ifdef PRINT
  printf("Sortiert.\n");
  fflush(stdout);
  #endif

  //Sukzessives Unterteilen:
  int *Unterteilung=Guarded_Malloc(sizeof(int)*signumber); //Unterteilung
  for(i=0;i<signumber;i++)Unterteilung[i]=0;
  int number=1;
  int number2=1;
  int drinne,draus;
  for(i=0;i<anzahl;i++) //Durch die C_Groups
  {
    if(Drop_Off[I[i]]<0.000001) // || i==0)
    {
      for(k=0;k<number;k++) // Alle bisherigen Aufteilungen werden auf Aufteilbarkeit untersucht:
      {
        drinne=0;
        draus=0;
        for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
        {
          if(Unterteilung[j]==k) //Ist in Aufteilung k
          {
            //printf("%d %d %d\n",i,I[i],j);
            //fflush(stdout);
            if(GrElement(C_Groups[I[i]],j)) {drinne+=1;}
            else {draus+=1;}
          }
        }
        //printf("gr %d -> drinne %d draus %d\n",Groupsize(C_Groups[I[i]]), drinne,draus);
        if(drinne>mingroup && draus>mingroup) //Wenn die Auteilung nicht mit bisherigen clasht wird aufgeteilt:
        {
          for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
          {
            if(Unterteilung[j]==k) //Ist in Aufteilung k
            {
              //printf("- %d %d %d\n",i,I[i],j);
              //fflush(stdout);
              if(GrElement(C_Groups[I[i]],j)) {Unterteilung[j]=number2;}
              else {Unterteilung[j]=number2+1;}
            }
          }
          number2+=2; 
        }
      }
      number=number2;
      number=Unterteilungskomprimierung(Unterteilung);
      Unterteilung_Assessment(Unterteilung);
    }
  }

  #ifdef PRINT
  printf("Sukzessives Unterteilen.\n");
  fflush(stdout);
  #endif

  //Jetzt kommt das relative Unterteilen:
  int c,count;
  int oldnumber=number+1;
  int bestgroup_i;
  double best_dropoff;
  unsigned long *Placeholder;
  int *Vars;

  int splitfactor;
  int * DieUnaufteilbaren=Guarded_Malloc(sizeof(int)*signumber);
  int unaufteilbar=0;
  int rep; //Der Repräsentant eines Clusters.

  while(oldnumber!=number)
  {
    oldnumber=number;
    for(k=0;k<number;k++)  //Durch die Aufteilungen
    {
      //Erstmal zähle ich, ob man die überhaupt aufteilen kann:
      count=0;
      for(i=0;i<signumber;i++){if(Unterteilung[i]==k){count++;rep=i;}}

      //Falls einer der Unaufteilbaren Reps zu diesem Cluster gehört, dann ist der Cluster unaufteilbar.
      for(i=0;i<unaufteilbar;i++)
      {
        if(Unterteilung[DieUnaufteilbaren[i]]==k)count=0;
      }

      if(count>mingroup*2)
      {
        //printf("count %d\n",count );
        best_dropoff=1000000.0;
        bestgroup_i=-1;

        for(i=0;i<anzahl;i++) //Durch die C_Groups 
        {
          c=Relative_Dropoff_Cutoff(I[i], 0, Unterteilung, k); //Relativer Drop_Off
          Placeholder=C_Groups[I[i]]; //Aufbewahren und später wiederreinschreiben.
          C_Groups[I[i]]=CliqueGroup(Cliques[I[i]], c);

          if(best_dropoff>Drop_Off[I[i]]) //Hier merke ich mir die beste Group falls sie unterteilt.
          {
            drinne=0;
            draus=0;
            for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
            {
              if(Unterteilung[j]==k) //Ist in Aufteilung k
              {
                //printf("%d %d %d\n",i,I[i],j);
                //fflush(stdout);
                if(GrElement(C_Groups[I[i]],j)) {drinne+=1;}
                else {draus+=1;}
              }
            }

            if(drinne>mingroup && draus>mingroup)
            {
              best_dropoff=Drop_Off[I[i]];
              bestgroup_i=i;
              //printf("best_dropoff %f\n",best_dropoff );
            }
          }

          //Unterteilen falls sinnvoll:
          if(Drop_Off[I[i]]<0.000001)
          {
            drinne=0;
            draus=0;
            for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
            {
              if(Unterteilung[j]==k) //Ist in Aufteilung k
              {
                //printf("%d %d %d\n",i,I[i],j);
                //fflush(stdout);
                if(GrElement(C_Groups[I[i]],j)) {drinne+=1;}
                else {draus+=1;}
              }
            }
            //printf("gr %d -> drinne %d draus %d\n",Groupsize(C_Groups[I[i]]), drinne,draus);
            if(drinne>mingroup && draus>mingroup) //Wenn die Auteilung nicht mit bisherigen clasht wird aufgeteilt:
            {
              for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
              {
                if(Unterteilung[j]==k) //Ist in Aufteilung k
                {
                  //printf("- %d %d %d\n",i,I[i],j);
                  //fflush(stdout);
                  if(GrElement(C_Groups[I[i]],j)) {Unterteilung[j]=number2;}
                  else {Unterteilung[j]=number2+1;}
                }
              }
              number2+=2; 
            }   
          }  
          free(C_Groups[I[i]]);
          C_Groups[I[i]]=Placeholder;
        }
        //Hier unterteile ich mit der besten Group:
        if(bestgroup_i>-1 && best_dropoff<0.000001)
        {
          drinne=0;
          draus=0;
          for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
          {
            if(Unterteilung[j]==k) //Ist in Aufteilung k
            {
              //printf("%d %d %d\n",i,I[i],j);
              //fflush(stdout);
              if(GrElement(C_Groups[I[bestgroup_i]],j)) {drinne+=1;}
              else {draus+=1;}
            }
          }
          //printf("gr %d -> drinne %d draus %d\n",Groupsize(C_Groups[I[i]]), drinne,draus);
          if(drinne>mingroup && draus>mingroup) //Wenn die Auteilung nicht mit bisherigen clasht wird aufgeteilt:
          {
            for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
            {
              if(Unterteilung[j]==k) //Ist in Aufteilung k
              {
                //printf("- %d %d %d\n",i,I[i],j);
                //fflush(stdout);
                if(GrElement(C_Groups[I[bestgroup_i]],j)) {Unterteilung[j]=number2;}
                else {Unterteilung[j]=number2+1;}
              }
            }
            number2+=2; 
          }   
        } 

/*        else if(0):  //Hier kommt die Kmeans mit rein:
        {
          Vars=Relative_Vars(Unterteilung, k, MaxCorrs, cutoff, mingroup);
          splitfactor=Kmeans(Unterteilung, k, Vars, mingroup);
          if(splitfactor==0)
          {
            DieUnaufteilbaren[unaufteilbar]=rep;
          }
          free(Vars);
        }*/
      }         
    }
    number=number2;
    number=Unterteilungskomprimierung(Unterteilung);
    Unterteilung_Assessment(Unterteilung);
  }

  //Es ist besser Kmeans separat laufen zu lassen:
  for(k=0;k<number;k++)  //Durch die Aufteilungen
  {
    count=0;
    for(i=0;i<signumber;i++){if(Unterteilung[i]==k){count++;}}
    if(count>mingroup*2)
    {
      Vars=Relative_Vars(Unterteilung, k, MaxCorrs, cutoff, mingroup);
      splitfactor=Kmeans(Unterteilung, k, Vars, mingroup);
      free(Vars);
    }
  }
  
  #ifdef PRINT
  printf("Unterteilt.\n");
  fflush(stdout);
  #endif

  return Unterteilung;
}


// Hier unterteile ich die Sigs mit relativen Dropoff-Cutoff-C_Groups und benutze nur noch Drop_Off==0
int *DropOff_Subdivision(double *MaxCorrs, double cutoff, double dropoffcutoff,int sizecutoff,int anfang, int ende, int mingroup)
{
  //Sortieren der C_Groups nach Drop_Off, Corr average und Size:
  int anzahl=0;
  int *I=Guarded_Malloc(sizeof(int)*siglength*5); //Indizes
  int i,j,k;
  for(i=0;i<siglength*5;i++)
  {
    if(MaxCorrs[i]>cutoff && Sizes[i]>sizecutoff)
    {
      I[anzahl]=i;
      anzahl++;
    }
  }

  #ifdef PRINT
  printf("Es gibt %d Groups über den Cutoffs.\n",anzahl );
  #endif

  for(i=0;i<anzahl;i++)
  {
    for(j=i+1;j<anzahl;j++)
    {
      if(Drop_Off[I[i]]>Drop_Off[I[j]]){k=I[i];I[i]=I[j];I[j]=k;}
      else if(Drop_Off[I[i]]==Drop_Off[I[j]])
      {
        if(Sizes[I[i]]<Sizes[I[j]]){k=I[i];I[i]=I[j];I[j]=k;}
        else if(Sizes[I[i]]==Sizes[I[j]])
        {
          if(MaxCorrs[I[i]]<MaxCorrs[I[j]]){k=I[i];I[i]=I[j];I[j]=k;}
        }
      }
    }
  }
  #ifdef PRINT
  printf("Sortiert.\n");
  fflush(stdout);
  #endif

  //Sukzessives Unterteilen:
  int *Unterteilung=Guarded_Malloc(sizeof(int)*signumber); //Unterteilung
  for(i=0;i<signumber;i++)Unterteilung[i]=0;
  int number=1;
  int number2=1;
  int drinne,draus;
  for(i=0;i<anzahl;i++) //Durch die C_Groups
  {
    if(Drop_Off[I[i]]<dropoffcutoff) // || i==0)
    {
      for(k=0;k<number;k++) // Alle bisherigen Aufteilungen werden auf Aufteilbarkeit untersucht:
      {
        drinne=0;
        draus=0;
        for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
        {
          if(Unterteilung[j]==k) //Ist in Aufteilung k
          {
            //printf("%d %d %d\n",i,I[i],j);
            //fflush(stdout);
            if(GrElement(C_Groups[I[i]],j)) {drinne+=1;}
            else {draus+=1;}
          }
        }
        //printf("gr %d -> drinne %d draus %d\n",Groupsize(C_Groups[I[i]]), drinne,draus);
        if(drinne>mingroup && draus>mingroup) //Wenn die Auteilung nicht mit bisherigen clasht wird aufgeteilt:
        {
          for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
          {
            if(Unterteilung[j]==k) //Ist in Aufteilung k
            {
              //printf("- %d %d %d\n",i,I[i],j);
              //fflush(stdout);
              if(GrElement(C_Groups[I[i]],j)) {Unterteilung[j]=number2;}
              else {Unterteilung[j]=number2+1;}
            }
          }
          number2+=2; 
        }
      }
      number=number2;
      number=Unterteilungskomprimierung(Unterteilung);
      Unterteilung_Assessment(Unterteilung);
    }
  }

  #ifdef PRINT
  printf("Sukzessives Unterteilen in %d.\n",number);
  fflush(stdout);
  #endif
  free(I);
  return Unterteilung;
}


void RelativeDropoff_Subdivision(int *Unterteilung,double *MaxCorrs, double cutoff, double dropoffcutoff,int sizecutoff,int anfang, int ende, int mingroup)
{
  int anzahl=0;
  int *I=Guarded_Malloc(sizeof(int)*siglength*5); //Indizes
  int i,j,k;
  for(i=0;i<siglength*5;i++)
  {
    if(MaxCorrs[i]>cutoff && Sizes[i]>sizecutoff)
    {
      I[anzahl]=i;
      anzahl++;
    }
  }

  int number=Unterteilungskomprimierung(Unterteilung);

  //Jetzt kommt das relative Unterteilen:
  int c,count;
  int bestgroup_i;
  double best_dropoff;
  unsigned long *Placeholder;
  int drinne,draus;

  for(k=0;k<number;k++)  //Durch die Aufteilungen
  {
    //Erstmal zähle ich, ob man die überhaupt aufteilen kann:
    count=0;
    for(i=0;i<signumber;i++)if(Unterteilung[i]==k){count++;}
    
    if(count>mingroup*2)
    {
      //printf("count %d\n",count );
      best_dropoff=1000000.0;
      bestgroup_i=-1;
      for(i=0;i<anzahl;i++) //Durch die C_Groups 
      {
        c=Relative_Dropoff_Cutoff(I[i], 0, Unterteilung, k); //Relativer Drop_Off
        Placeholder=C_Groups[I[i]]; //Aufbewahren und später wiederreinschreiben.
        C_Groups[I[i]]=CliqueGroup(Cliques[I[i]], c);
        if(best_dropoff>Drop_Off[I[i]]) //Hier merke ich mir die beste Group falls sie unterteilt.
        {
          drinne=0;
          draus=0;
          for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
          {
            if(Unterteilung[j]==k) //Ist in Aufteilung k
            {
              //printf("%d %d %d\n",i,I[i],j);
              //fflush(stdout);
              if(GrElement(C_Groups[I[i]],j)) {drinne+=1;}
              else {draus+=1;}
            }
          }
          if(drinne>mingroup && draus>mingroup)
          {
            best_dropoff=Drop_Off[I[i]];
            bestgroup_i=i;
            //printf("best_dropoff %f\n",best_dropoff );
          }
        }
        
        //Unterteilen falls sinnvoll:
        if(Drop_Off[I[i]]<dropoffcutoff)
        {
          drinne=0;
          draus=0;
          for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
          {
            if(Unterteilung[j]==k) //Ist in Aufteilung k
            {
              //printf("%d %d %d\n",i,I[i],j);
              //fflush(stdout);
              if(GrElement(C_Groups[I[i]],j)) {drinne+=1;}
              else {draus+=1;}
            }
          }
          //printf("gr %d -> drinne %d draus %d\n",Groupsize(C_Groups[I[i]]), drinne,draus);
          if(drinne>mingroup && draus>mingroup) //Wenn die Auteilung nicht mit bisherigen clasht wird aufgeteilt:
          {
            for(j=0;j<signumber;j++) //Die werden jetzt differenziert:
            {
              if(Unterteilung[j]==k) //Ist in Aufteilung k
              {
                //printf("- %d %d %d\n",i,I[i],j);
                //fflush(stdout);
                if(GrElement(C_Groups[I[i]],j)) {Unterteilung[j]=number+1+k*2;}
                else {Unterteilung[j]=number+2+k*2;}
              }
            }
          }   
        }  
        free(C_Groups[I[i]]);
        C_Groups[I[i]]=Placeholder;
      }
    }         
  }

  number=Unterteilungskomprimierung(Unterteilung);
  Unterteilung_Assessment(Unterteilung);
    
  #ifdef PRINT
  printf("Rel. Unterteilt in %d.\n",number);
  fflush(stdout);
  #endif 
}


// Hier unterteile ich die Sigs mit relativen Dropoff-Cutoff-C_Groups und benutze nur noch Drop_Off==0
void Kmeans_Subdivision(int *Unterteilung,double *MaxCorrs, double cutoff, int mingroup)
{
  int *Vars;
  int number=Unterteilungskomprimierung(Unterteilung);
  int i,k,splitfactor,count;
  for(k=0;k<number;k++)  //Durch die Aufteilungen
  {
    count=0;
    for(i=0;i<signumber;i++){if(Unterteilung[i]==k){count++;}}
    if(count>mingroup*2)
    {
      Vars=Relative_Vars(Unterteilung, k, MaxCorrs, cutoff, mingroup);
      splitfactor=Kmeans(Unterteilung, k, Vars, mingroup);
      free(Vars);
    }
  }
  number=Unterteilungskomprimierung(Unterteilung);    
  #ifdef PRINT
  printf("Unterteilt in %d.\n",number);
  fflush(stdout);
  #endif
}


void Results_und_Connections(int *Unterteilung)
{
  int i;
  //
  for(i=0;i<signumber;i++)printf("%d -> %d\n",i,Unterteilung[i]);
  // Hier wird vielleicht der Rest unterteilt?:
  Unterteilung_Assessment(Unterteilung);
  
  int *Unterteilung1=Guarded_Malloc(sizeof(int)*signumber);
  int *Unterteilung2=Guarded_Malloc(sizeof(int)*signumber);
  for(i=0;i<signumber;i++)
  {
    if(i%30<15){Unterteilung1[i]=i/30;Unterteilung2[i]=-1;}
    else {Unterteilung2[i]=i/30;Unterteilung1[i]=-1;}
  }
  Verbindungen(Unterteilung1,Unterteilung,Unterteilung2);

  int **Unterteilungen=Guarded_Malloc(sizeof(int)*10);
  int unt_no=3;
  int ii;

  //Multi_Verbindungen(Unterteilungen,unt_no);
  for(unt_no=3;unt_no<10;unt_no++)
  {
    Unterteilungen[0]=Unterteilung1;
    Unterteilungen[1]=Unterteilung;
    Unterteilungen[unt_no-1]=Unterteilung2;

    for(i=2;i<unt_no-1;i++)
    {
      Unterteilungen[i]=Guarded_Malloc(sizeof(int)*signumber);
      for(ii=0;ii<signumber;ii++)Unterteilungen[i][ii]=Unterteilung[(ii+30*(i-1))%signumber];
    }

    //Multi_Verbindungen(Unterteilungen,unt_no);
    printf("\n\nBack_And_Forth(Unterteilungen,unt_no);: %d\n",unt_no);
    Back_And_Forth(Unterteilungen,unt_no);
  }
}

















int A_C(int i,int j, int t)
{
  int schnitt=Schnitt(Groups[i],Groups[j]); 
  int gr1=Schnitt(Groups[i],LocalCoverage[j/5]); 
  int gr2=Schnitt(LocalCoverage[i/5],Groups[j]); 
  if(gr1<gr2)
  {
    if(schnitt+t>gr1)
    {
      return 1;
    }
  }
  return 0;
}

// Hier analysiere ich mal herum: 
void A_C_Ana(int i,int j, int k)
{
  int l;
  int histo[100];
  for(l=0;l<100;l++)histo[l]=0;

  printf("Gr1 Elemente:\n");
  for(l=0;l<signumber;l++)
  {
    if(GrElement(Groups[i],l))
    {
      printf("%d ",l );
      histo[l/30]++;
    }
  }
  printf("\n");

  printf("Copygruppen: ");
  for(l=0;l<100;l++)if(histo[l]>5)printf("%d-%d=%d ",l*30,(l+1)*30,histo[l] );
  printf("\n");

  printf("Gr2^Gr3 nicht Gr1\n");
  for(l=0;l<signumber;l++)
  {
    if(GrElement(Groups[j],l) && GrElement(Groups[k],l) && !GrElement(Groups[i],l))printf("%d ",l );
  }
  printf("\n");

  printf("Gr1 nicht Gr3^Gr2\n");
  for(l=0;l<signumber;l++)
  {
    if(GrElement(Groups[i],l) && !(GrElement(Groups[k],l) || GrElement(Groups[j],l)))printf("%d ",l );
  }
  printf("\n");

  printf("\n");
}

// Hier schaue ich mir mal approximate containment an
void Approximate_Containment_Structure(double *MaxCorrs, unsigned long **Groups, double cutoff, int mincov)
{
  int i,j,k;

  for(i=10000;i<siglength*5;i++)
  {
    if(MaxCorrs[i]>cutoff && Groupsizearray[i]<100 && Coverage[i/5]==3000)
    {
      for(j=10000;j<siglength*5;j++)
      { 
        if(MaxCorrs[j]>cutoff && Groupsizearray[j]<100)
        {
          if(Schnitt(LocalCoverage[i/5],LocalCoverage[j/5])==Groupsize(LocalCoverage[i/5]) && A_C(i,j,mincov))
          {
            for(k=10000;k<siglength*5;k++)
            {
              if(MaxCorrs[k]>cutoff && Groupsizearray[k]<300)
              {
                if(Schnitt(LocalCoverage[i/5],LocalCoverage[k/5])==Groupsize(LocalCoverage[k/5]) )
                {
                  if(Schnitt(LocalCoverage[j/5],LocalCoverage[k/5])==Groupsize(LocalCoverage[j/5]) )
                  {
                    if(A_C(j,k,mincov))
                    {
                      if(!A_C(i,k,mincov))
                      {
                        if(Schnitt(Groups[i],Groups[k])*3>Groupsizearray[i])
                        {
                          printf("%d < %d < %d\n",Groupsizearray[i],Groupsizearray[j],Groupsizearray[k]);
                          A_C_Ana(i,j,k);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}











//Bis hier sind die Funktionen jetzt ohne Sig.

void Clique_Aussieber(double *MaxCorrs, double cutoff, int min_cliquesize, int min_schnitt, int sig)
{
  int i,j;
  int count=0;

  //Hier sieben wir mal Cliquen aus, die eventuell inkorrekt sind. 
  UnderCut=Guarded_Malloc(sizeof(int)*siglength*5);
  for(i=0;i<siglength*5;i++)UnderCut[i]=0;
  for(i=0;i<siglength*5;i++)
  {
    if(MaxCorrs[i]>cutoff && Sizes[i]>min_cliquesize)
    {
      count++;
      for(j=i+1;j<siglength*5;j++)
      {
        if(MaxCorrs[j]>cutoff && Sizes[j]>min_cliquesize)
        {
          if(Schnitt(C_Groups[i], C_Groups[j])<min_schnitt)
          {
            UnderCut[i]++;
            UnderCut[j]++;
          }
        }
      }
    }
    else
    {
      MaxCorrs[i]=0.0; //Besser schon mal rauskicken.
    }
  } 
  #ifdef PRINT
  printf("Hier habe ich %d potentielle Cliquen.\n",count );
  #endif

  //Jetzt entfernen wir nach und nach die schlechtesten bis alle gleich gut sind:
  int *ausgesiebt=Guarded_Malloc(sizeof(int)*siglength*5);
  count=0;
  for(i=0;i<siglength*5;i++)
  {
    if((MaxCorrs[i]>cutoff && Sizes[i]>min_cliquesize)){ausgesiebt[i]=0;count++;}
    else {ausgesiebt[i]=1;}
  }
  #ifdef PRINT
  printf("Ja, tatsächlich %d potentielle Cliquen.\n",count );
  #endif

  int bestundercut=10000000;
  int worstundercut=0;
  int worst_i=0;
  count=2;

  while(bestundercut!=worstundercut && count>1)
  {
    bestundercut=1000000;
    worstundercut=0;    
    worst_i=-1;
    //bestes und schlechtestes berechnen
    for(i=0;i<siglength*5;i++)
    {
      if(ausgesiebt[i]==0)
      {
        if(UnderCut[i]<bestundercut){bestundercut=UnderCut[i];}
        if(UnderCut[i]>worstundercut){worstundercut=UnderCut[i];worst_i=i;}
      }
    }
    //schlechtestes aussieben
    if(bestundercut<worstundercut){ausgesiebt[worst_i]=1;}
    #ifdef PRINT
    printf("worst_i: %d, %d < %d ?\n",worst_i, bestundercut,worstundercut);
    #endif
    //neues UnderCut berechnen
    count=0;
    if(worst_i>-1)
    {
      for(j=0;j<siglength*5;j++)
      {
        if(ausgesiebt[j]==0){if(Schnitt(C_Groups[worst_i], C_Groups[j])<min_schnitt){UnderCut[j]--;}count++;}
      }
      #ifdef PRINT
      printf("Count %d\n",count );
      #endif
    }
  }

  //Kleiner Test:
  for(i=0;i<siglength*5;i++)UnderCut[i]=0;
  for(i=0;i<siglength*5;i++)
  {
    if(ausgesiebt[i]==0)
    {
      for(j=i+1;j<siglength*5;j++)
      {
        if(ausgesiebt[j]==0)
        {
          if(Schnitt(C_Groups[i], C_Groups[j])<min_schnitt)
          {
            UnderCut[i]++;
            UnderCut[j]++;
          }
        }
      }
    }
  } 
  #ifdef PRINT
  for(i=0;i<siglength*5;i++)
  {
    if(UnderCut[i]>0)printf("%d ",UnderCut[i]);
  }
  printf("\n");
  #endif

  //Wir manipulieren MaxCorrs als output
  count=0;
  for(i=0;i<siglength*5;i++)
  {
    if(ausgesiebt[i]!=0)MaxCorrs[i]=0.0;
    else count++;
  }
  #ifdef PRINT
  printf("Es gibt noch %d Cliquen\n",count);
  #endif
  
}

// Hier verfeinere ich die letzten Groups durch Schnitte:
void Clique_Schnitte(double *MaxCorrs, double cutoff, int sig)
{
  int i,j,jj,positives,falsepositives;
  int min_undercut=1000000;
  for(i=0;i<siglength*5;i++){if(MaxCorrs[i]>cutoff){if(UnderCut[i]<min_undercut){min_undercut=UnderCut[i];}}}
  
  int minsize=1000000;
  j=0;
  for(i=0;i<siglength*5;i++)
  {
    if(MaxCorrs[i]>cutoff)
    {
      if(minsize>Groupsize(C_Groups[i]) && UnderCut[i]==min_undercut )
      {
        j=i;
        minsize=Groupsize(C_Groups[i]);
      }
    }
  }

  unsigned long *SchnittGroup=NULL;
  SchnittGroup=GrInitialize();
  GrCopy(C_Groups[j],SchnittGroup);
  minsize=Groupsize(SchnittGroup);

  MaxCorrs[j]=0.0; //Die erste Gruppe wird ausgesiebt.

  #ifdef PRINT
  printf("Jetzt die Schneiderei...\n");
  fflush(stdout);
  #endif

  int schnitt;
  for(j=0;j<5;j++)
  {
    jj=-1;
    for(i=0;i<siglength*5;i++)
    {
      if(MaxCorrs[i]>cutoff)
      {
        schnitt=Schnitt(SchnittGroup,C_Groups[i]);
        if(schnitt<minsize && schnitt>10)
        {
          minsize=schnitt;
          jj=i;
        }
      }
    }

    if(jj>-1)
    {
      for(i=0;i<signumber;i++)
      {
        if(!GrElement(C_Groups[jj],i))
        {
          GrDel(SchnittGroup,i);
        }
      }
      MaxCorrs[jj]=0.0; //Die Gruppe wird ausgesiebt.
    }
    else
    {
      j=1000;
    }
  }
  //Hier überprüfe ich wie gut der Schnitt ist:
  falsepositives=0;
  positives=0;
  for(i=0;i<signumber;i++)
  {
    if(GrElement(SchnittGroup,i))
    {
      if(i>=(sig/30)*30 && i<((sig/30)+1)*30)
      {positives++;}
      else 
      {falsepositives++;}
    }
  }
  #ifdef PRINT
  printf("Schnitt %d %d\n",minsize,Groupsize(SchnittGroup));
  printf("%d positives, %d falsepositives.\n",positives,falsepositives); 
  #endif

}


// Hier gehe ich durch die un-signifikanten Groups und berechne ihre Fehlerrate:
void ErrorEstimate(double *MaxCorrs, double cutoff)
{
  int i,j,jj;
  int grsize[5];
  int cov;
  double basebasecount=0.0;
  double basebaseerror=0.0;
  double basegapcount=0.0;
  double basegaperror=0.0;
  double gapbasecount=0.0;
  double gapbaseerror=0.0;

  for(i=0;i<siglength;i++)
  {
    if(MaxCorrs[i*5+0]<cutoff && MaxCorrs[i*5+1]<cutoff && MaxCorrs[i*5+2]<cutoff && MaxCorrs[i*5+3]<cutoff && MaxCorrs[i*5+4]<cutoff) //Column ohne Variation
    {
      for(j=0;j<5;j++){grsize[j]=Groupsize(Groups[i*5+j]);printf("%d ",grsize[j] );}
      printf("\n");
      cov=grsize[0]+grsize[1]+grsize[2]+grsize[3]+grsize[4];
      for(j=0;j<4;j++) //Basengroups
      {
        if(grsize[j]>intmax(intmax(grsize[(j+1)%5],grsize[(j+2)%5]),intmax(grsize[(j+3)%5],grsize[(j+4)%5])))
        {
          for(jj=0;jj<4;jj++)
          {
            if(jj!=j) //Die kleinen Basengruppen sind Error:
            {
              basebaseerror+=(double)grsize[jj]/(double)cov;
            }
          }
          basebasecount+=3.0;
          //Die Gapgruppe ist auch error
          basegaperror+=(double)grsize[4]/(double)cov;
          basegapcount+=1.0;
        }
      }
      if(grsize[4]>intmax(intmax(grsize[0],grsize[1]),intmax(grsize[2],grsize[3]))) //Gapgroups
      {
        for(jj=0;jj<4;jj++)
        {
          gapbaseerror+=(double)grsize[jj]/(double)cov;
        } 
        gapbasecount+=4.0;       
      }
    }
  }
  printf("In Basencolumns sind %f falsche Basen und %f falsche Gaps.\n",basebaseerror/basebasecount,basegaperror/basegapcount );
  printf("In Gapcolumns sind %f falsche Basen.\n",gapbaseerror/gapbasecount);
}








void Help()
{
  printf("Usage: ./RepeatResolver MSA_path\n");
  printf("Flags:\n");
  printf("-f start end    Start and end column of the resolution window.\n");
  printf("-c <x>          Expected sequencing coverage.\n");
  printf("-t <x>          The threshold for the significance of base group intersections.\n");
  exit(0);
}












int main(int argc, char *argv[])
{ 

sc_p=&sc;
char *MApath_p;

if(argc<2){printf("Usage: ./RepeatResolver MApath <options>\n");exit(0);}

MApath_p=argv[1];
int cov=30;
int i;
char *covstring_p;
int parallel=0;
int von=-1;
int bis=-1;
double cutoff=0.0;

if(argc>1)  // Für weitere options
{
  char *output_p;
  for(i=2;i<argc;i++)
  {
    if(argv[i][0]=='-' && argv[i][1]=='o')
    {
      if(i+1<argc)output_p=argv[i+1];
    }  
  }
  for(i=2;i<argc;i++)
  {
    if(argv[i][0]=='-' && argv[i][1]=='p')
    {
      if(i+1<argc)
      {
        parallel=strtol(argv[i+1],(char **)NULL,10);
        printf("NTHREADS: %d\n",parallel);
      }
    }  
  }
  for(i=2;i<argc;i++)
  {  
    if(argv[i][0]=='-' && argv[i][1]=='h')
    {
      Help();
    } 
  }
  for(i=2;i<argc;i++)
  {

    if(argv[i][0]=='-' && argv[i][1]=='c')
    {
      if(i+1<argc)
      {
        covstring_p=argv[i+1];
        cov=strtol(covstring_p,(char **)NULL,10);
        printf("Coverage %d\n",cov );
      }
    }  
  }  

  for(i=2;i<argc;i++)
  {

    if(argv[i][0]=='-' && argv[i][1]=='f')
    {
      if(i+2<argc)
      {
        von=strtol(argv[i+1],(char **)NULL,10);
        bis=strtol(argv[i+2],(char **)NULL,10);
        printf("Full coverage from column %d until %d.\n",von, bis);
      }
    }  
  }  
  
  for(i=2;i<argc;i++)
  {
    if(argv[i][0]=='-' && argv[i][1]=='t')  //t for threshold
    {
      if(i+1<argc)
      {
      cutoff=atof(argv[i+1]);
      }
    }  
  } 
}

if(von==-1 && bis==-1)
{
  von=0;
  bis=Max_Var_Anzahl;
}

Einlesen(MApath_p,von,bis);
//exit(0);


// Hier suchen wir einen generischen Path für die MaxCorrs und berechnen sie wenn er nicht existiert.
// Hier berechne ich die MaxCorrs, wenn sie nicht schon existieren:
char CorrNames[300];

char von_string[10];
sprintf(von_string, "%d", von);
char bis_string[10];
sprintf(bis_string, "%d", bis);

strcpy(CorrNames,"MaxCorrsOf_");
strcat(CorrNames,MApath_p);

printf("%s\n",CorrNames);

double *MaxCorrs=MaxCorrsEinlesen(CorrNames,von,bis); 


printf("Opened \n");fflush(stdout);

if(cutoff<0.1)cutoff=-1.0*log10(1.0/(double)(siglength*5.0));


int significant=0;
for(i=0;i<siglength*5;i++){if(MaxCorrs[i]>cutoff)significant++;}
printf("%d correlations make the cutoff.\n",significant );
//exit(0);
printf("Cutoff %f\n", cutoff);fflush(stdout);


if(MaxCorrs==NULL)
{
  printf("AllMaxCorrsRechner\n");
  if(parallel)  //parallelized
  {
    MaxCorrs=Parallel_AllMaxCorrsRechner(parallel, 0, siglength*5, cov, signumber, cutoff);
  }
  else //Single thread
  {
    MaxCorrs=AllMaxCorrsRechner(0, siglength*5, cov, signumber, cutoff);
  }
  MaxCorrsRausschreiben(MaxCorrs, CorrNames);
}



// Coverageeinschränkung: Vielleicht später ein bisschen robuster
int maxcov=0;
for(i=0;i<siglength;i++)
{
  if(Coverage[i]>maxcov)maxcov=Coverage[i];
}
printf("Maxcov: %d\n",maxcov );

for(i=0;i<siglength*5;i++)
{
  if(Coverage[i/5]*10<maxcov*9)MaxCorrs[i]=0.0;
}

// Dann Group_Refinement
double greedy=cutoff;
int anfang=0;  // column 0
int ende=siglength; // letzte column
int mincov=cov;
int maxclique=30;

if(parallel){Parallel_Group_Refinement(MaxCorrs, cutoff, anfang, ende, mincov, maxclique, greedy,parallel);}
else{Group_Refinement(MaxCorrs, cutoff, anfang, ende, mincov, maxclique, greedy);}

// Und die Perfect Subdivision:
int sizecutoff=-1;
int mingroup=mincov/2;

/*time_t start_time=time(NULL);
int *Unterteilung=Perfect_Subdivision(MaxCorrs, cutoff, sizecutoff, anfang, ende, mingroup);
time_t end_time=time(NULL);
printf("Runtime: %lu sec.\n", (end_time-start_time));*/

time_t start_time=time(NULL);
double dropoffcutoff=0.0001;
int *Unterteilung=DropOff_Subdivision(MaxCorrs, cutoff, dropoffcutoff, sizecutoff, anfang, ende, mingroup);
time_t end_time=time(NULL);
printf("DropOff_Subdivision Runtime: %lu sec.\n", (end_time-start_time));
char UnterteilungNames[300];
strcpy(UnterteilungNames,"DropoffSubdivisionOf_");
strcat(UnterteilungNames,von_string);
strcat(UnterteilungNames,"_");
strcat(UnterteilungNames,bis_string);
strcat(UnterteilungNames,"_");
strcat(UnterteilungNames,MApath_p);
int *Komp_Unterteilung=UnterteilungsKomplettierung(Unterteilung);
Unterteilung_Rausschreiben(Komp_Unterteilung, UnterteilungNames);


start_time=time(NULL);
RelativeDropoff_Subdivision(Unterteilung,MaxCorrs, cutoff, dropoffcutoff, sizecutoff, anfang, ende, mingroup);
end_time=time(NULL);
printf("RelativeDropoff_Subdivision Runtime: %lu sec.\n", (end_time-start_time));
strcpy(UnterteilungNames,"RelDropSubdivisionOf_");
strcat(UnterteilungNames,von_string);
strcat(UnterteilungNames,"_");
strcat(UnterteilungNames,bis_string);
strcat(UnterteilungNames,"_");
strcat(UnterteilungNames,MApath_p);
Komp_Unterteilung=UnterteilungsKomplettierung(Unterteilung);
Unterteilung_Rausschreiben(Komp_Unterteilung, UnterteilungNames);

start_time=time(NULL);
Kmeans_Subdivision(Unterteilung, MaxCorrs, cutoff, mingroup);
end_time=time(NULL);
printf("KmeansSubdivision Runtime: %lu sec.\n", (end_time-start_time));
strcpy(UnterteilungNames,"KmeansSubdivisionOf_");
strcat(UnterteilungNames,von_string);
strcat(UnterteilungNames,"_");
strcat(UnterteilungNames,bis_string);
strcat(UnterteilungNames,"_");
strcat(UnterteilungNames,MApath_p);
Komp_Unterteilung=UnterteilungsKomplettierung(Unterteilung);
Unterteilung_Rausschreiben(Komp_Unterteilung, UnterteilungNames);

exit(0);
Results_und_Connections(Unterteilung);

exit(0);



}



//gcc -Wall -I/usr/local/include -c RepeatResolver.c
//gcc -L/usr/local/lib RepeatResolver.o -lgsl -lgslcblas -lm -o RepeatResolver -lpthread
//./RepeatResolver MMA_path -c 30 -f von bis





