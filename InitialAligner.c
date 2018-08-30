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

#include <pthread.h>


//#define PRINT


/*
This tool creates an initial multiple sequence alignment by aligning reads to a template.
The resulting msa is only a rough arrangement of the most important bases and has to be refined.

This version is parallelized to manage very large MSAs. 
*/


void *Guarded_Malloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL)
    { fprintf(stderr,"\nError in InitialAligner\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}

void *Guarded_Realloc(void *p, size_t size)
{ p = realloc(p,size);
  if (p == NULL)
    { fprintf(stderr,"\nError in InitialAligner:\n");
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


/* Hier ein neuer Plan wie ich das ganze parallelisiere:

- Einmal durch die Reads gehen und offsets speichern.
- ReadingFasta liest spezifischen Read in spezifisches Memory
- Hilfsfunktion alloziert Memory für die Matrix
- Und geht dann durch alle x-te Reads und füllt den Alignment-Array
- Parallelfunktion steckt Hilfsfunktion in verschiedene threads

*/

int **Alignments;
int *Readlengths;
double *AlignmentError;
off_t *Offsets;
int ReadCounter(char *Path)
{
  char *buffer;
  buffer=Guarded_Malloc(sizeof(char)*70000);
  char *s;
  FILE * File=NULL;
  size_t len=70000;
  File=fopen(Path,"r");
  int count=0;
  if(File==NULL){exit(1);}
  while((s = fgets(buffer, len, File)) != NULL)
  {
    if(buffer[0]=='>')count++;
  }
  fclose(File);
  Alignments=Guarded_Malloc(sizeof(int*)*count);
  Readlengths=Guarded_Malloc(sizeof(int)*count);
  AlignmentError=Guarded_Malloc(sizeof(double)*count);
  Offsets=Guarded_Malloc(sizeof(off_t)*count);
  free(buffer);
  return count;
}


int Offsetter(char *Path)
{

  char *buffer;
  buffer=Guarded_Malloc(sizeof(char)*70000);

  char *s;
  FILE * File=NULL;
  size_t len=70000;

  File=fopen(Path,"r");
  int i;
  int readcount=0;
  off_t offset=0;

  while(!feof(File))
  {
    s = fgets(buffer, len, File);

    if(s==NULL)
    {
      if(feof(File)){}
      else{printf("Error in ReadCounter\n");exit(1);}
    }

    else if(buffer[0]=='>')
    {
      Offsets[readcount]=offset;
      readcount++;
      i=0;
      while(buffer[i]!='\n')  //The beginning of the read
      {
        i+=1;
      }
      offset+=i+1;  //+1 wegen \n
    }
    else  //The bases are read in 
    {
      //printf("\n\n%d\n",readcount);
      i=0;
      while(buffer[i]!='\n')
      {
        i++;
/*        if(buffer[i]=='A' || buffer[i]=='a')i++;
        else if(buffer[i]=='C' || buffer[i]=='c')i++;
        else if(buffer[i]=='G' || buffer[i]=='g')i++;
        else if(buffer[i]=='T' || buffer[i]=='t')i++;*/
      }
      offset+=i+1;  //+1 wegen \n
    }

  }

  fclose(File);
  free(buffer);


  return readcount;

}

int ReadingFasta(char *Path, int readcount, char *Read, char *buffer)
{
  char *s;
  FILE * File=NULL;
  size_t len=70000;

  File=fopen(Path,"r");
  if(fseek(File,Offsets[readcount],SEEK_SET)){printf("The End.\n");fflush(stdout);exit(1);}

  int i;

  //printf("%p\n",File );
  fflush(stdout);

  if(File==NULL){printf("File == NULL.\n");fflush(stdout);exit(1);}

  int basecount=0;
  int bigger=0;
  while(bigger<2)
  {
    s = fgets(buffer, len, File);
    //printf("#%c#%d#\n",buffer[0],readcount );
    if(s==NULL)
    {
      if(feof(File)){bigger+=1;}
      else{printf("Error in ReadingFasta.\n");fflush(stdout); exit(0);}
    }

    else if(buffer[0]=='>')
    {
      bigger+=1;
      if(bigger==2)
      {
        fclose(File);
        Read[basecount]='\0';
        return basecount;        
      }
    }
    else  //The bases are read in 
    {
      i=0;
      while(buffer[i]!='\n')
      {
        
        if(buffer[i]=='A' || buffer[i]=='a')Read[basecount]='a';
        else if(buffer[i]=='C' || buffer[i]=='c')Read[basecount]='c';
        else if(buffer[i]=='G' || buffer[i]=='g')Read[basecount]='g';
        else if(buffer[i]=='T' || buffer[i]=='t')Read[basecount]='t';
        else {basecount--;}       
        basecount++;
        i++;
      }
    }
  }

  fclose(File);
  Read[basecount]='\0';
  return basecount;

}



//Reading in the template
char Template[70000];
int templatelength=0;
void ReadingTemplate(char *Path)
{
  char *buffer=Guarded_Malloc(sizeof(char)*70000);
  char *s;
  FILE * File=NULL;
  size_t len=70000;
  File=fopen(Path,"r");

  int i;

  if(File==NULL){exit(1);}

  int basecount=0;

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
  free(buffer);
  Template[basecount]='\n';
  templatelength=basecount;

}






void Help()
{
  printf("Usage: ./InitialAligner_parallel template.fasta Seq.fasta\n");
  printf("Flags:\n");
  printf("-o msa_path    Path of the resulting multiple sequence alignment. Default: SimulatedMSA.\n");
  printf("-s <150>       Path of the seq class information: Which seqs are instances of the repeat. Default: SimulatedSeqClass\n");
  printf("-e <0.30>      The mapping error cutoff being used to detect instances of the template.\n");
  exit(0);
}

//The aligner with matrix and co as argument. 
int IntoAligner(char *shortstring, char *longstring, int length1, int length2, char **Matrix, unsigned long *Row, int aligncount,char *EditScript)
{
  #ifdef PRINT
  printf("Into IntoAligner %d %d \n",length1,length2);fflush(stdout);
  #endif
  //int length1=strlen(shortstring);
  //int length2=strlen(longstring);
  int x,y;
  unsigned long UpperRow;
  unsigned long NewEntry;
  for(y=-1;y<length2;y++)Row[y]=0;

  //Matrix filling: 0 is diagonal=y-1 and x-1; 1 is up = y-1,x; 2 is left=y,x-1
  int m;
  for(x=0;x<length1;x++)
  {
    UpperRow=x;
    Row[-1]=x+1;
    for(y=0;y<length2;y++)
    {

      m=1;
      if(shortstring[x]==longstring[y])m=0;

      NewEntry=UpperRow+m;
      if(m){Matrix[x][y]=0;}
      else{Matrix[x][y]=3;}

      if(Row[y-1]+1<NewEntry)  
      {
        NewEntry=Row[y-1]+1;
        Matrix[x][y]=1;
      }

      if(Row[y]+1<NewEntry)
      {
        NewEntry=Row[y]+1;
        Matrix[x][y]=2;
      }
      UpperRow=Row[y];
      Row[y]=NewEntry;
    }
  }

  #ifdef PRINT
  printf("Matrix done\n");fflush(stdout);
  #endif

  // Backtracking: 
  int count=0;

  x=length1-1;
  y=length2-1;

  int min,i,einstieg_y;
  einstieg_y=y;
  min=Row[einstieg_y];
  for(i=length2-1;i>0;i--)
  {
    if(Row[i]<min)
    {
      min=Row[i];
      einstieg_y=i;
    }
  }

  #ifdef PRINT
  printf("Entry %d\n",einstieg_y );
  printf("%lu\n",Row[einstieg_y]);
  #endif

  AlignmentError[aligncount]=(double)Row[einstieg_y]/(double)length1;

  #ifdef PRINT
  printf("Alignment Error[%d]=%f\n",aligncount,AlignmentError[aligncount] );
  #endif

  while(y>einstieg_y)
  {
    EditScript[count]='i';
    count++;
    y--;
  }

  while(x>-1 && y>-1)
  {
    //printf("%d %d\n",x,y );
    //printf("%d %d %d -> %d\n",Matrix[x-1][y-1],Matrix[x][y-1],Matrix[x-1][y],Matrix[x][y] );
    m=1;  
    if(shortstring[x]==longstring[y])m=0;

    if(Matrix[x][y]==0){EditScript[count]='s';count++;x--;y--;}
    else if(Matrix[x][y]==3){EditScript[count]='m';count++;x--;y--;}
    else if(Matrix[x][y]==1){EditScript[count]='i';count++;y--;}
    else if(Matrix[x][y]==2){EditScript[count]='d';count++;x--;}
    else
    {
      printf("Backtracking error.\n");
      exit(1);
    } 
    //printf("x:%d y:%d Matr:%d\n",x,y,Matrix[x][y]);
  }

  #ifdef PRINT
  printf("Backtracking done\n");fflush(stdout);
  #endif  

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

  //Inverting edit script:
  for(x=0;x<count/2;x++)
  {
    m=EditScript[x];
    EditScript[x]=EditScript[count-x-1];
    EditScript[count-x-1]=m;
  }

  //for(x=0;x<count;x++)printf("%c",EditScript[x]);
  //printf("\n"); 

  // Turning the edit script into a alignment:

  Alignments[aligncount]=Guarded_Malloc(sizeof(int)*length1);
  Readlengths[aligncount]=length1;

  #ifdef PRINT
  printf("Readlengths[%d]=%d\n",aligncount,Readlengths[aligncount] );
  printf("EditScript length=%d\n",count );
  #endif

  x=0;
  y=0;
  for(i=0;i<count;i++)
  {
    if(EditScript[i]=='s' || EditScript[i]=='m')
    {
      //printf("%c - %c\n", shortstring[x],longstring[y]);
      Alignments[aligncount][x]=y;
      x++;
      y++;
    }
    if(EditScript[i]=='i')
    {
      //printf("_ - %c\n", longstring[y]);    
      y++;  
    }
    if(EditScript[i]=='d')
    {
      //printf("%c - _\n", shortstring[x]);     
      Alignments[aligncount][x]=-1;
      x++;
    }
  }
  //for(x=0;x<length1;x++)printf("%d ",Alignments[aligncount][x] );
  //printf("\n");
  #ifdef PRINT
  printf("Out IntoAligner\n");fflush(stdout);
  #endif

  return min;
}

char *Readspath_p;
//This function aligns every n-th read to the template:
void *All_Aligner(void *x)
{
  int * args;
  args=((int *) x);

  int start=(*(args+0));
  int step=(*(args+1)); 
  int maxlength1=(*(args+2));
  int maxlength2=(*(args+3));
  int readanzahl=(*(args+4));

  //Now we allocate an alignment matrix for this set of reads:
  unsigned long *RowSchatten=Guarded_Malloc(sizeof(unsigned long)*(maxlength2+1));
  unsigned long *Row=&RowSchatten[1];
  char **Matrix=Guarded_Malloc(sizeof(char*)*maxlength1);
  char *EditScript=Guarded_Malloc(sizeof(char)*(maxlength1+maxlength2));
  char *Read=Guarded_Malloc(sizeof(char)*maxlength1);
  char *buffer=Guarded_Malloc(sizeof(char)*maxlength1*2);
  int i;
  for(i=0;i<maxlength1;i++)Matrix[i]=Guarded_Malloc(sizeof(char)*maxlength2);
  int readlength;
  // Here we step through the Reads
  int aligncount;
  int percentage=1;
  for(aligncount=start;aligncount<readanzahl;aligncount+=step)
  {
    if(start==0)
    {
      if((aligncount*100)/readanzahl>percentage){printf("%d percent.\n",percentage );percentage+=1;}
    }
    //Reading in the read
    //printf("Into ReadingFasta\n");fflush(stdout);
    readlength=ReadingFasta(Readspath_p, aligncount, Read, buffer);
    //printf("Outo ReadingFasta\n");fflush(stdout);
    //Aligning
    IntoAligner(Read, Template, readlength, maxlength2, Matrix, Row, aligncount, EditScript);
  }
  return NULL;
}

// Hier ist Code von dem erfolgreich parallelisiertem Repeatresolver:



void Parallel_Aligning(int maxlength1, int maxlength2, int NTHREADS, int readanzahl)
{

  pthread_t threads[NTHREADS];
  int *thread_args[NTHREADS];

  int rc,i;

  #ifdef PRINT
  time_t start_time=time(NULL);
  printf("Parallel_Aligning\n");
  #endif


  // spawn the threads 
  for (i=0; i<NTHREADS; ++i)
  {
    //(int anfang, int ende, int mincov, int maxgroup, double cutoff, int NTHREADS, int thread)
    thread_args[i] = Guarded_Malloc(sizeof(int)*5);
    thread_args[i][0]=i;
    thread_args[i][1]=NTHREADS;
    thread_args[i][2]=maxlength1;  
    thread_args[i][3]=maxlength2;
    thread_args[i][4]=readanzahl;

    #ifdef PRINT
    printf("spawning thread %d\n", i);
    #endif
    rc = pthread_create(&threads[i], NULL, All_Aligner, (void *) thread_args[i]);
  }

  #ifdef PRINT
  printf("threads done\n");fflush(stdout);
  #endif

  // wait for threads to finish 
  for (i=0; i<NTHREADS; ++i) {
    rc = pthread_join(threads[i], NULL);
  }
  #ifdef PRINT
  printf("pthread_join\n");fflush(stdout);
  #endif

  #ifdef PRINT
  time_t time_spent=time(NULL)-start_time;
  printf("%ld sec.\n", time_spent);  
  #endif

}



void Building_MSA(char *outputfile, char *seqclassfile, char *Readspath_p, double errorcutoff, int readanzahl, int templatelength, int maxlength1)
{
  FILE *datei;
  datei=fopen(outputfile,"w");

  FILE *datei2;
  datei2=fopen(seqclassfile,"w");

  //printf("Files opened\n");fflush(stdout);

  char *Read=Guarded_Malloc(sizeof(char)*maxlength1);
  char *buffer=Guarded_Malloc(sizeof(char)*maxlength1*2);
  int readlength;

  //Alignment[j][i] tells us where the i-th base of the j-th read went. 
  // ==-1 indicates that the i-th base is placed between two template bases.

  //Determine the maximal gap between all bases:
  int *Gapcount=Guarded_Malloc(sizeof(int)*templatelength+1);
  int i,j,k,l,count,gap;
  for(i=0;i<templatelength+1;i++)Gapcount[i]=0;

  for(j=0;j<readanzahl;j++)
  {
    gap=0;
    count=0;
    //first gap 
    i=0;
    while(Alignments[j][i]==-1)i++;
    gap=Alignments[j][i];  //The bases before 'gap'

    for(i=0;i<Readlengths[j];i++)
    {
      if(Alignments[j][i]==-1)
      {
        count++;
        if(count>Gapcount[gap]){Gapcount[gap]=count;}
      }
      else
      {
        gap=Alignments[j][i]+1;  //The bases before Alignments[j][i] are already counted.
        count=0;
      }
    }
  }

  //printf("Gaps determined\n");fflush(stdout);

  //Write out the reads while respecting the gaps:
  for(j=0;j<readanzahl;j++)
  {
    //printf("(");fflush(stdout);
    readlength=ReadingFasta(Readspath_p, j, Read, buffer);
    //printf("%d,%d) ",j,readlength );fflush(stdout);
    if(AlignmentError[j]<errorcutoff)
    {
      fprintf(datei2,"r\n");
      if(readlength>0)
      {
        k=0;
        i=0;
        while(i<templatelength+1)
        {
          //gap bases into the gap
          count=0;
          while(k<readlength && Alignments[j][k]==-1)
          {
            fprintf(datei,"%c",Read[k]);
            if(Read[k]=='\n'){printf("Bumm\n"); exit(0);}
            k++;
            count++;
          }
          //Filling the gap with '-'
          for(l=count;l<Gapcount[i];l++)fprintf(datei,"-");

          //Next base base if it exists
          if(k<readlength && Alignments[j][k]==i)
          {
            fprintf(datei,"%c",Read[k]);
            if(Read[k]=='\n'){printf("Bumm\n"); exit(0);}
            k++;
          }
          else
          {
            fprintf(datei,"-");
          }
          i++;
        }
      }
      else
      {
        for(i=0;i<templatelength+1;i++)
        {
          fprintf(datei,"-");
          for(k=0;k<Gapcount[i];k++)fprintf(datei,"-");
        }
      }
      fprintf(datei,"\n"); 
    }
    //if(j==410)for(i=0;i<readlength;i++)printf("%d\n", Alignments[j][i]);

    else //No repeat seq, 'l' is for large scale variation as opposed to 'r' for repetitive. 
    {
      fprintf(datei2,"l\n");
    }
  }
  printf("\nFiles written.\n");fflush(stdout);

  fclose(datei);
  fclose(datei2);
}




int main(int argc, char *argv[])
{
  if(argc<3)Help();

  Readspath_p=argv[2];
  char *Templatepath_p=argv[1];


  //Here the outputpath is specified 
  //Building a filename with the right prefix Type+'_05perc_30kb_Template.fasta
  char outputfile[300];
  char seqclassfile[300];

  int i=0;
  char TemplateSuffix[]="Template.fasta";
  char outputsuffix[]="MSA";
  char seqclasssuffix[]="SeqClass";
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

  strcpy(seqclassfile,DataPrefix);
  strcat(seqclassfile,seqclasssuffix);

  //////////////////////////
  char *seqclass_p=&seqclassfile[0];
  char *output_p=&outputfile[0];

  double errorcutoff=0.30;
  int NTHREADS=1;



  for(i=1;i<argc;i++)
  {

    if(argv[i][0]=='-' && argv[i][1]=='o')
    {
      //printf("%s\n",argv[i]);
      if(i+1<argc)output_p=argv[i+1];
    }  

    if(argv[i][0]=='-' && argv[i][1]=='s')
    {
      //printf("%s\n",argv[i]);
      if(i+1<argc)seqclass_p=argv[i+1];
    }  

    if(argv[i][0]=='-' && argv[i][1]=='e')
    {
      errorcutoff=atof(argv[i+1]);
      printf("errorcutoff %f.\n",errorcutoff );
    }   

    if(argv[i][0]=='-' && argv[i][1]=='p')
    {
      NTHREADS=atoi(argv[i+1]);
    } 

    if(argv[i][0]=='-' && argv[i][1]=='h')
    {
      Help();
    } 
  }

  if(argc<2){printf("Usage: ./InitialAligner template.fasta Seq.fasta\n");exit(0);}



  ReadingTemplate(Templatepath_p);
  printf("template length %d\n",templatelength );

  int maxlength2=templatelength;
  int maxlength1=40000;  //maximal readlength


  printf("output file: %s\n",output_p );
  printf("seqclass file: %s\n",seqclass_p );

  int readanzahl=ReadCounter(Readspath_p);  //Also allocates memory

  printf("read count %d\n",readanzahl );

  Offsetter(Readspath_p);

  Parallel_Aligning(maxlength1, maxlength2, NTHREADS, readanzahl);

  printf("Writing the msa.\n");
  Building_MSA(output_p, seqclass_p, Readspath_p, errorcutoff, readanzahl, templatelength, maxlength1);
  exit(0);
}  



