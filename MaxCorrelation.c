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

#define Max_Var_Anzahl 150000
#define Max_Sig_Anzahl 30000

#define PRINT



/*
Hier berechne ich nur die MaxCorrs für ein MSA. 
*/



unsigned long **Groups;
int *Groupsizearray;
unsigned long **LocalCoverage;

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
int *Coverage;

void Einlesen(char* MApath_p, int von, int bis)   // aka Input füllen.
{ 
  char *s;  
  int i,j,k;
  char buffer[Max_Var_Anzahl];

 
  //Statt char Signatures[Max_Sig_Anzahl][Max_Var_Anzahl];  mal mallocen
  Signatures=Guarded_Malloc(sizeof(char*)*Max_Sig_Anzahl);

  signumber=-1;

  FILE * File;
  File=fopen(MApath_p,"r");
  if(File==NULL){printf("MA is missing.\n"); exit(1);}


  while ((s=fgets(buffer, Max_Var_Anzahl-2, File)) != NULL)
  {
    if(signumber==-1)
    {
      siglength = strlen(buffer)-1;
      //Now we know the siglength, i.e. the number of variations and can allocate memory:
      Groups=Guarded_Malloc(sizeof(unsigned long*)*siglength*5);
      Groupsizearray=Guarded_Malloc(sizeof(unsigned long*)*siglength*5);
      LocalCoverage=Guarded_Malloc(sizeof(unsigned long*)*siglength);
      Coverage=Guarded_Malloc(sizeof(int)*siglength);
    }

    if(strlen(buffer)-1==siglength)   //buffer[von]!=' ' && buffer[bis]!=' ')  //We try to calculate all Corr first look at sections with full cov later
    {
      signumber++;
      *(Signatures+signumber)=Guarded_Malloc(sizeof(char)*siglength);

      for(i=0;i<siglength;i++)
      {
        if(buffer[i]=='a' || buffer[i]=='A')
        {
          Signatures[signumber][i]=0;
        }
        else if(buffer[i]=='c' || buffer[i]=='C')
        {
          Signatures[signumber][i]=1;
        }
        else if(buffer[i]=='g' || buffer[i]=='G')
        {
          Signatures[signumber][i]=2;
        }
        else if(buffer[i]=='t' || buffer[i]=='T')
        {
          Signatures[signumber][i]=3;
        }
        else if(buffer[i]=='-' || buffer[i]=='_')
        {
          Signatures[signumber][i]=4;
        }
        else
        {
          Signatures[signumber][i]=5;
        }
      }      
    }
  }


  signumber++;  //Damit es nicht der letzte Index sondern die Anzahl ist.
  printf("There are %d sequences.\n",signumber);
  printf("Siglength is %d.\n",siglength);

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
  //for(ii=anfang;ii<ende;ii++)
  for(ii=0;ii<siglength;ii++)
  {
    //printf("%d\n",ii);fflush(stdout);
    for(k=0;k<5;k++)
    {
      i=ii*5+k; //Die relevante Var in Column ii
      //printf("%d\n", Groupsizearray[i]);
      if(Groupsizearray[i]>mincov/4 && Groupsizearray[i]<maxgroup)
      {
        //for(jj=ii+20;jj<ende;jj++) 
        for(jj=ii+20;jj<siglength;jj++)  
        {
          //Wenn der Abstand zu groß wird gibt es keine gemeinsame coverage mehr, dann wird abgebrochen:
          //printf("%d %d %d\n",ii,jj,siglength);fflush(stdout);
          if(Schnitt(LocalCoverage[ii],LocalCoverage[jj])<mincov)
          {
            //printf("%d\n",Schnitt(LocalCoverage[ii],LocalCoverage[jj]) );
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
                printf("%f\n",Z );
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

double *HilfsMaxCorr[128];
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
  int baseno;

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
  //for(ii=116170;ii<ende;ii++)
  {
    if(thread==0)
    {
      //printf("%d ",ii );fflush(stdout);
      if( (percentage*(ende-anfang))/100 == ii )
      {
        current_time=time(NULL);
        #ifdef PRINT
        //printf("%d %% in %ld min -> %ld min.\n",percentage,(current_time-start_time)/60,(((current_time-start_time)*(100-percentage))/(percentage*60))  );  
        printf("%d %% (%d/%d) in %ld min of %ld min\n",percentage,ii,ende,(current_time-start_time)/60, (10000*(current_time-start_time)/ ( (200-percentage)*percentage ))/60);
        #endif
        percentage++;
      }
      fflush(stdout); 
    }     


    if(ii%NTHREADS==thread)
    {
      baseno=Groupsizearray[ii*5]+Groupsizearray[ii*5+1]+Groupsizearray[ii*5+2]+Groupsizearray[ii*5+3];
      for(k=0;k<5;k++)
      {
        i=ii*5+k; //Die relevante Var in Column ii
        if(Groupsizearray[i]>mincov/4 && Groupsizearray[i]<maxgroup && baseno>Coverage[ii]/2)
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
                  //printf("%d %d %d %d %d In \n",j,i,jj,ii,ende);fflush(stdout);
                  Z=PositiveSignificance(i,j);
                  //printf("%d %d %d %d %d Out \n",j,i,jj,ii,ende) ;fflush(stdout);
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

  int rc,i,j;

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
  printf("threads spawned\n");fflush(stdout);
  #endif

  /* wait for threads to finish */
  for (i=0; i<NTHREADS; ++i) {
    rc = pthread_join(threads[i], NULL);
  }

  #ifdef PRINT
  printf("pthread_join %d\n",rc);fflush(stdout);
  #endif

  // MaxCorrs vereinigen:
  for(i=0;i<siglength*5;i++)
  {
    for(j=1;j<NTHREADS;j++)
    {
      if(HilfsMaxCorr[j][i]>HilfsMaxCorr[0][i])
      {
        HilfsMaxCorr[0][i]=HilfsMaxCorr[j][i];
      }
    }
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







int main(int argc, char *argv[])
{ 

sc_p=&sc;
char *MApath_p;

if(argc<2){printf("Usage: ./MaxCorrelation MSApath <options>\n");exit(0);}

MApath_p=argv[1];
int cov=30;
int i;
char *covstring_p;
int parallel=1;
int von=-1;
int bis=-1;

if(argc>1)  // Für weitere options
{

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
}

time_t start_time=time(NULL);

Einlesen(MApath_p,von,bis);

//exit(0);
if(von==-1 && bis==-1)
{
  von=0;
  bis=siglength;
}
printf("From %d to %d\n",von,bis );


// Hier suchen wir einen generischen Path für die MaxCorrs und berechnen sie wenn er nicht existiert.
// Hier berechne ich die MaxCorrs, wenn sie nicht schon existieren:
char CorrNames[300];
strcpy(CorrNames,"MaxCorrsOf_");
strcat(CorrNames,MApath_p);

printf("%s\n",CorrNames );
double *MaxCorrs=NULL;  //=MaxCorrsEinlesen(CorrNames,von,bis); 

double cutoff=-1.0*log10(1.0/(double)(siglength*5.0));

printf("Cutoff %f\n", cutoff);fflush(stdout);


if(MaxCorrs==NULL)
{
  printf("AllMaxCorrs\n");
  if(parallel)  //parallelized
  {
    MaxCorrs=Parallel_AllMaxCorrsRechner(parallel, 0, siglength, cov, signumber, cutoff);
  }
  else //Single thread
  {
    MaxCorrs=AllMaxCorrsRechner(0, siglength, cov, signumber, cutoff);
  }
  MaxCorrsRausschreiben(MaxCorrs, CorrNames);
}

time_t end_time=time(NULL);
printf("Runtime: %lu sec.\n", (end_time-start_time));

exit(0);  //Erstmal nur die MaxCorrs berechnen.





}



//gcc -Wall -I/usr/local/include -c MaxCorrelation.c
//gcc -L/usr/local/lib MaxCorrelation.o -lgsl -lgslcblas -lm -o MaxCorrelation -lpthread
//./MaxCorrelation MMA_path -c 30 -p 4





