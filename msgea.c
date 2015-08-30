#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#define GSLH 0

#define NSTAT 7
#define NSPEC 3
#define NGO 2100000
#define ALPHA 0.01
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

typedef struct {
	char exist;
	char temp;
	int level;
	double pvalue[NSTAT];
//	unsigned long tempcount[NSPEC];
//	unsigned long tempbackground[NSPEC];
	unsigned long count[NSPEC];
	unsigned long background[NSPEC];
	unsigned long icount[NSPEC];
	unsigned long ibackground[NSPEC];
	unsigned int n_parents;
	unsigned int n_childs;
	unsigned long *parents;
	unsigned long *childs;
} go_term;

typedef struct {
	unsigned long count[NSPEC];
	unsigned long background[NSPEC];
} go_counts;


go_term graph[NGO];
unsigned long n_genes_background[NSPEC],n_genes[NSPEC];

#if GSLH == 1
double my_cdf_hypergeometric_Q (unsigned int k, unsigned int n1, unsigned int n2, unsigned int t){
  return gsl_cdf_hypergeometric_Q(k,n1,n2,t);
}
double my_ran_hypergeometric_pdf (unsigned int k, unsigned int n1, unsigned int n2, unsigned int t){
  return gsl_ran_hypergeometric_pdf(k,n1,n2,t);
}
#else
double my_cdf_hypergeometric_Q (unsigned int k, unsigned int n1, unsigned int n2, unsigned int t){
  return gsl_cdf_binomial_Q(k, (double)t/(double)(n1+n2), n1);
  //gsl_cdf_hypergeometric_Q(k,n1,n2,t);
}
double my_ran_hypergeometric_pdf (unsigned int k, unsigned int n1, unsigned int n2, unsigned int t){
  return gsl_ran_binomial_pdf(k, (double)t/(double)(n1+n2), n1);
  //gsl_ran_hypergeometric_pdf(k,n1,n2,t);
}
#endif

double gsl_cdf_hypergeometric_t (unsigned int k, unsigned int n1, unsigned int n2, unsigned int t){
  return my_cdf_hypergeometric_Q(k,n1,n2,t)+my_ran_hypergeometric_pdf(k,n1,n2,t);
}


int check_getline(ssize_t c){
  if ((c>0) && (c>0)) {
    return 1;
  } else {
    return -1;
  };
};

void init_graph(){
  unsigned long j,k;
  unsigned int s,g_line;
  ssize_t bytes_read;
  char *line,*split;
  size_t n_line;
  
  FILE * file_graph;
//  graph=malloc(NGO*sizeof(struct go_tem));
  for(j=1;j<=NGO;j++){
    graph[j-1].exist=0;
    graph[j-1].temp=0;
    graph[j-1].level=-1;
    graph[j-1].n_parents=0;
    graph[j-1].n_childs=0;
    for(s=1;s<=NSPEC;s++){
//      graph[j-1].tempcount[s-1]=0;
//      graph[j-1].tempbackground[s-1]=0;
      graph[j-1].count[s-1]=0;
      graph[j-1].background[s-1]=0;
      graph[j-1].icount[s-1]=0;
      graph[j-1].ibackground[s-1]=0;
    };
    for(s=1;s<=NSTAT;s++){
      graph[j-1].pvalue[s-1]=1;
    };
  };
  
  line=malloc(2*sizeof(char));
  n_line=2;
  
  file_graph=fopen("msgea_dag.dat","r");
  for(;check_getline(bytes_read=getline(&line,&n_line,file_graph))>0;){
    g_line=(unsigned int)(bytes_read/11);
    split=strtok(line,"\t");
    sscanf(split,"GO:%lu",&j);
    graph[j-1].exist=1;
    for(s=2;s<=g_line;s++){
      split=strtok(NULL,"\t");
      sscanf(split,"GO:%lu",&k);
      graph[j-1].n_parents++;
      graph[k-1].n_childs++;
      if(graph[j-1].n_parents==1){
	graph[j-1].parents=(unsigned long *)malloc(sizeof(unsigned long));
      } else {
	graph[j-1].parents=(unsigned long *)realloc(graph[j-1].parents,sizeof(unsigned long)*graph[j-1].n_parents);
      };
      graph[j-1].parents[graph[j-1].n_parents-1]=k;
      if(graph[k-1].n_childs==1){
	graph[k-1].childs=(unsigned long *)malloc(sizeof(unsigned long));
      } else {
	graph[k-1].childs=(unsigned long *)realloc(graph[k-1].childs,sizeof(unsigned long)*graph[k-1].n_childs);
      };
      graph[k-1].childs[graph[k-1].n_childs-1]=j;
      };
  }; 
  fclose(file_graph);

  
};

int get_level(unsigned long j){
  int l,i,t;
  if (graph[j-1].n_parents==0){
    return 1;
  } else {
    l=1000;
    for (i=1;i<=graph[j-1].n_parents;i++){
      if ((t=1+get_level(graph[j-1].parents[i-1]))<l) { l=t; };
    };
    return l;
  };
};

void init_levels(){
  unsigned int j;
  for(j=1;j<=NGO;j++){
    if (graph[j-1].exist==1) { graph[j-1].level=get_level(j); };
  };
};

void init_background(char *filename, int species){
  FILE * file_bg;
  unsigned long j,k;
  unsigned int s,g_line;
  ssize_t bytes_read;
  char *line,*split;
  size_t n_line;

  line=malloc(2*sizeof(char));
  n_line=2;
  
  n_genes_background[species-1]=0;
  file_bg=fopen(filename,"r");
  for(;check_getline(bytes_read=getline(&line,&n_line,file_bg))>0;){
    sscanf(line,"GO:%lu",&k);
    graph[k-1].background[species-1]++;
    graph[k-1].ibackground[species-1]++;
    n_genes_background[species-1]++;
  };
  fclose(file_bg);
};


void init_counts(char *filename, int species){
  FILE * file_bg;
  unsigned long j,k;
  unsigned int s,g_line;
  ssize_t bytes_read;
  char *line,*split;
  size_t n_line;

  line=malloc(2*sizeof(char));
  n_line=2;
  
  n_genes[species-1]=0;
  file_bg=fopen(filename,"r");
  if (file_bg == NULL) {
        fprintf(stderr, "Can't open input file %s. Program terminated.\n", filename);
        exit(2);
  }
  for(;check_getline(bytes_read=getline(&line,&n_line,file_bg))>0;){
    sscanf(line,"GO:%lu",&k);
    graph[k-1].count[species-1]++;
    graph[k-1].icount[species-1]++;
    n_genes[species-1]++;
  };
  fclose(file_bg);
};


go_counts get_child_count(unsigned long child){
  unsigned int s,c;
  go_counts cb,cb2;
  for(s=1;s<=NSPEC;s++){
   cb.count[s-1]=0;
   cb.background[s-1]=0;
  };
  if((graph[child-1].exist==0)||(graph[child-1].temp==1)) {

    return cb;
  };
  graph[child-1].temp=1;
  for(s=1;s<=NSPEC;s++){
    cb.count[s-1]+=graph[child-1].icount[s-1];
    cb.background[s-1]+=graph[child-1].ibackground[s-1];
  }; 
  if (graph[child-1].n_childs>0) { 
    for(c=1;c<=graph[child-1].n_childs;c++){
      cb2=get_child_count(graph[child-1].childs[c-1]);
      for(s=1;s<=NSPEC;s++){
	cb.count[s-1]+=cb2.count[s-1];
	cb.background[s-1]+=cb2.background[s-1];
      };
    };
  };
  return cb;
};

void reset_temp(){
  unsigned long i;
  unsigned int s;
  for(i=1;i<=NGO;i++){
//    for(s=1;s<=NSPEC;s++){
     if (graph[i-1].exist==1){ graph[i-1].temp=0; };
//      graph[i-1].tempcount[s-1]=0;
//      graph[i-1].tempbackground[s-1]=0;
//    };
  };
};

void reset_temp2(unsigned long j){
  unsigned long i;
  unsigned int c;
  if(j>0){
    graph[j-1].temp=0; 
    for(c=1;c<=graph[j-1].n_childs;c++){
      reset_temp2(graph[j-1].childs[c-1]);
    };
  };
};


void propagate(){
  unsigned long i,j;
  unsigned int s,c;
  go_counts cb;
  j=0;
  for(i=1;i<=NGO;i++){
    if (graph[i-1].exist==1){
    reset_temp2(j);
    j=i;
    for(c=1;c<=graph[i-1].n_childs;c++){
      cb=get_child_count(graph[i-1].childs[c-1]);
      for(s=1;s<=NSPEC;s++){
	graph[i-1].count[s-1]+=cb.count[s-1];
	graph[i-1].background[s-1]+=cb.background[s-1];
      };
    };
    };
  };
};


void compute_stat(){
  double p, stat, stat_temp, stat_temp2, den[NSPEC];
  unsigned long i;
  unsigned int s, sm, k1, k2, k3;
  unsigned int kam,kbm,kcm;

  for(i=1;i<=NGO;i++){
    if (graph[i-1].exist==1){
    //printf("node %u\n",i);
    printf(".");
    sm=0;
    for (s=1;s<=NSPEC;s++){
      if (graph[i-1].count[s-1]>0) {
	if (graph[i-1].background[s-1]>=graph[i-1].count[s-1]) {
	graph[i-1].pvalue[s-1]=gsl_cdf_hypergeometric_t(graph[i-1].count[s-1], n_genes[s-1], n_genes_background[s-1]-n_genes[s-1], graph[i-1].background[s-1]);      
	sm++;
	} else {graph[i-1].pvalue[s-1]=2;};
      };
    };
    //printf("Single species statistics\n");
    if ((sm>0)&&(i!=3674)&&(i!=5575)&&(i!=8150)){
    sm=0;
    p=0;
    stat=0;
    stat_temp2=0;
    for (s=1;s<=NSPEC;s++){
      if (graph[i-1].background[s-1]>0) { 
	den[s-1]=sqrt((double)((n_genes[s-1])*(graph[i-1].background[s-1]))/(double)(n_genes_background[s-1])); 
      } else {
	if (graph[i-1].count[s-1]>0) {printf("something very strange in node %lu species %u\n",i,s); };
	den[s-1]=1;
      };
      stat_temp=(double)(graph[i-1].count[s-1])/den[s-1];
      stat+=stat_temp;
      if (stat_temp>stat_temp2) {sm=s; stat_temp2=stat_temp;};
    };
    if(stat<=0){
      for (s=1;s<=NSPEC*(NSPEC-1)/2+1;s++){
	graph[i-1].pvalue[NSPEC+s-1]=3;
      };
    } else if ((den[1-1]!=1)&&(den[2-1]!=1)&&(den[3-1]!=1)) {
      
      double p12,p13,p23,stat12,stat13,stat23;
      p12=0;
      p13=0;
      p23=0;
      stat12=(double)(graph[i-1].count[1-1])/den[1-1]+(double)(graph[i-1].count[2-1])/den[2-1];
      stat13=(double)(graph[i-1].count[1-1])/den[1-1]+(double)(graph[i-1].count[3-1])/den[3-1];
      stat23=(double)(graph[i-1].count[2-1])/den[2-1]+(double)(graph[i-1].count[3-1])/den[3-1];
      
      s=1;
      kam=(unsigned int)ceil(stat*den[s-1]);
      if (kam<0) kam=0;
      s=2;
      kbm=(unsigned int)ceil(stat*den[s-1]);
      if (kbm<0) kbm=0;
      s=3;
      kcm=(unsigned int)ceil(stat*den[s-1]);
      if (kcm<0) kcm=0;

      p+=gsl_cdf_hypergeometric_t(kam, n_genes[1-1], n_genes_background[1-1]-n_genes[1-1], graph[i-1].background[1-1])+
         gsl_cdf_hypergeometric_t(kbm, n_genes[2-1], n_genes_background[2-1]-n_genes[2-1], graph[i-1].background[2-1])+
	 gsl_cdf_hypergeometric_t(kcm, n_genes[3-1], n_genes_background[3-1]-n_genes[3-1], graph[i-1].background[3-1])+
       -(gsl_cdf_hypergeometric_t(kam, n_genes[1-1], n_genes_background[1-1]-n_genes[1-1], graph[i-1].background[1-1])*
         gsl_cdf_hypergeometric_t(kbm, n_genes[2-1], n_genes_background[2-1]-n_genes[2-1], graph[i-1].background[2-1])+
	 gsl_cdf_hypergeometric_t(kam, n_genes[1-1], n_genes_background[1-1]-n_genes[1-1], graph[i-1].background[1-1])*
         gsl_cdf_hypergeometric_t(kcm, n_genes[3-1], n_genes_background[3-1]-n_genes[3-1], graph[i-1].background[3-1])+
	 gsl_cdf_hypergeometric_t(kbm, n_genes[2-1], n_genes_background[2-1]-n_genes[2-1], graph[i-1].background[2-1])*
	 gsl_cdf_hypergeometric_t(kcm, n_genes[3-1], n_genes_background[3-1]-n_genes[3-1], graph[i-1].background[3-1]))
      +  gsl_cdf_hypergeometric_t(kam, n_genes[1-1], n_genes_background[1-1]-n_genes[1-1], graph[i-1].background[1-1])*
         gsl_cdf_hypergeometric_t(kbm, n_genes[2-1], n_genes_background[2-1]-n_genes[2-1], graph[i-1].background[2-1])*
	 gsl_cdf_hypergeometric_t(kcm, n_genes[3-1], n_genes_background[3-1]-n_genes[3-1], graph[i-1].background[3-1]);      
      for(k3=0;(k3+3<=kam+kbm+kcm)&&(p<ALPHA);k3++){
      for(k2=0;(k2+2<=min(k3+2,kam+kbm))&&(p<ALPHA);k2++){
      for(k1=0;(k1+1<=min(k2+1,kam))&&(p<ALPHA);k1++){
	if ((double)k1/den[0]+(double)(k2-k1)/den[1]+(double)(k3-k2)/den[2]>stat) p+=my_ran_hypergeometric_pdf(k1,n_genes[0],n_genes_background[0]-n_genes[0], graph[i-1].background[0])*
	  my_ran_hypergeometric_pdf(k2-k1,n_genes[1],n_genes_background[1]-n_genes[1], graph[i-1].background[1])*
	  my_ran_hypergeometric_pdf(k3-k2,n_genes[2],n_genes_background[2]-n_genes[2], graph[i-1].background[2]);
      };
      };
      };
      
      s=1;
      kam=(unsigned int)ceil(stat12*den[s-1]);
      if (kam<0) kam=0;
      s=2;
      kbm=(unsigned int)ceil(stat12*den[s-1]);
      if (kbm<0) kbm=0;
      s=3;
      kcm=(unsigned int)ceil(stat12*den[s-1]);
      if (kcm<0) kcm=0;
	 
      p12+=gsl_cdf_hypergeometric_t(kam, n_genes[1-1], n_genes_background[1-1]-n_genes[1-1], graph[i-1].background[1-1])+
         gsl_cdf_hypergeometric_t(kbm, n_genes[2-1], n_genes_background[2-1]-n_genes[2-1], graph[i-1].background[2-1])+
       -(gsl_cdf_hypergeometric_t(kam, n_genes[1-1], n_genes_background[1-1]-n_genes[1-1], graph[i-1].background[1-1])*
         gsl_cdf_hypergeometric_t(kbm, n_genes[2-1], n_genes_background[2-1]-n_genes[2-1], graph[i-1].background[2-1]));      
      for(k2=0;(k2+2<=kam+kbm)&&(p12<ALPHA);k2++){
      for(k1=0;(k1+1<=min(k2+1,kam))&&(p12<ALPHA);k1++){
	if ((double)k1/den[0]+(double)(k2-k1)/den[1]>stat12) p12+=my_ran_hypergeometric_pdf(k1,n_genes[0],n_genes_background[0]-n_genes[0], graph[i-1].background[0])*
	  my_ran_hypergeometric_pdf(k2-k1,n_genes[1],n_genes_background[1]-n_genes[1], graph[i-1].background[1]);
      };
      };

      s=1;
      kam=(unsigned int)ceil(stat13*den[s-1]);
      if (kam<0) kam=0;
      s=2;
      kbm=(unsigned int)ceil(stat13*den[s-1]);
      if (kbm<0) kbm=0;
      s=3;
      kcm=(unsigned int)ceil(stat13*den[s-1]);
      if (kcm<0) kcm=0;
      
      p13+=gsl_cdf_hypergeometric_t(kam, n_genes[1-1], n_genes_background[1-1]-n_genes[1-1], graph[i-1].background[1-1])+
         gsl_cdf_hypergeometric_t(kcm, n_genes[3-1], n_genes_background[3-1]-n_genes[3-1], graph[i-1].background[3-1])+
       -(gsl_cdf_hypergeometric_t(kam, n_genes[1-1], n_genes_background[1-1]-n_genes[1-1], graph[i-1].background[1-1])*
         gsl_cdf_hypergeometric_t(kcm, n_genes[3-1], n_genes_background[3-1]-n_genes[3-1], graph[i-1].background[3-1]));      
      for(k2=0;(k2+2<=kam+kcm)&&(p13<ALPHA);k2++){
      for(k1=0;(k1+1<=min(k2+1,kam))&&(p13<ALPHA);k1++){
	if ((double)k1/den[0]+(double)(k2-k1)/den[2]>stat13) p13+=my_ran_hypergeometric_pdf(k1,n_genes[0],n_genes_background[0]-n_genes[0], graph[i-1].background[0])*
	  my_ran_hypergeometric_pdf(k2-k1,n_genes[2],n_genes_background[2]-n_genes[2], graph[i-1].background[2]);
      };
      };
      
      s=1;
      kam=(unsigned int)ceil(stat23*den[s-1]);
      if (kam<0) kam=0;
      s=2;
      kbm=(unsigned int)ceil(stat23*den[s-1]);
      if (kbm<0) kbm=0;
      s=3;
      kcm=(unsigned int)ceil(stat23*den[s-1]);
      if (kcm<0) kcm=0;
      
      p23+=gsl_cdf_hypergeometric_t(kbm, n_genes[2-1], n_genes_background[2-1]-n_genes[2-1], graph[i-1].background[2-1])+
         gsl_cdf_hypergeometric_t(kcm, n_genes[3-1], n_genes_background[3-1]-n_genes[3-1], graph[i-1].background[3-1])+
       -(gsl_cdf_hypergeometric_t(kbm, n_genes[2-1], n_genes_background[2-1]-n_genes[2-1], graph[i-1].background[2-1])*
         gsl_cdf_hypergeometric_t(kcm, n_genes[3-1], n_genes_background[3-1]-n_genes[3-1], graph[i-1].background[3-1]));      
      for(k2=0;(k2+2<=kbm+kcm)&&(p23<ALPHA);k2++){
      for(k1=0;(k1+1<=min(k2+1,kbm))&&(p23<ALPHA);k1++){
	if ((double)k1/den[1]+(double)(k2-k1)/den[2]>stat23) p23+=my_ran_hypergeometric_pdf(k1,n_genes[1],n_genes_background[1]-n_genes[1], graph[i-1].background[1])*
	  my_ran_hypergeometric_pdf(k2-k1,n_genes[2],n_genes_background[2]-n_genes[2], graph[i-1].background[2]);
      };
      };

       
      
//      p+= gsl_cdf_hypergeometric_t(graph[i-1].count[s-1], n_genes[s-1], n_genes_background[s-1]-n_genes[s-1], graph[i-1].background[s-1]) ;      

      
      
      
      
      
      
      
/*      
      //rewrite
      for(k3=0;(k3<=graph[i-1].background[0]+graph[i-1].background[1]+graph[i-1].background[2])&&(k3<=n_genes[0]+n_genes[1]+n_genes[2]);k3++){
      for(k2=0;k2<=min(k3,min(graph[i-1].background[0]+graph[i-1].background[1],n_genes[0]+n_genes[1]));k2++){
      for(k1=0;k1<=min(k2,min(graph[i-1].background[0],n_genes[0]));k1++){
	printf("par %u %u %u\n",k1,k2,k3);
      //for(k2=0;(k2<=graph[i-1].background[0]+graph[i-1].background[1])&&(k2<=n_genes[0]+n_genes[1]);k2++){
      //for(k1=0;(k1<=graph[i-1].background[0])&&(k1<=n_genes[0]);k1++){
        if ((double)k1/den[0]+(double)(k2-k1)/den[1]+(double)(k3-k2)/den[2]>stat) p+=my_ran_hypergeometric_pdf(k1,n_genes[0],n_genes_background[0]-graph[i-1].background[0], graph[i-1].background[0])*
	  my_ran_hypergeometric_pdf(k2-k1,n_genes[1],n_genes_background[1]-graph[i-1].background[1], graph[i-1].background[1])*
	  my_ran_hypergeometric_pdf(k3-k2,n_genes[2],n_genes_background[2]-graph[i-1].background[2], graph[i-1].background[2]);
      if (p>ALPHA) break;
      };
      if (p>ALPHA) break;
      };
      if (p>ALPHA) break;
      };
      graph[i-1].pvalue[NSPEC]=p;
    */
      
      
      
      graph[i-1].pvalue[3]=p12;
      graph[i-1].pvalue[4]=p13;
      graph[i-1].pvalue[5]=p23;
      graph[i-1].pvalue[6]=p;
      
      
    } else {graph[i-1].pvalue[3]=2; graph[i-1].pvalue[4]=2; graph[i-1].pvalue[5]=2; graph[i-1].pvalue[6]=2;  };
    
    };
    
    };  
  };
};


void print_stat(char *filename, int is_3){
  unsigned long i,j;
  FILE * file_out;
  file_out=fopen(filename,"w");
  fprintf(file_out,"## Results of Multi-Species Gene Enrichment Analysis: \n ## Number of set and background GOs for each species: ");
  for (j=1;j<=NSPEC;j++){
    fprintf(file_out,"\t%lu\t%lu",n_genes[j-1],n_genes_background[j-1]);
  };
  fprintf(file_out,"\n");
  if(is_3==1){
      fprintf(file_out,"Go_term\tLevel\tOccurrences_species1\tBackground_species1\tOccurrences_species2\tBackground_species2\tOccurrences_species3\tBackground__species3\tPvalue_species1\tPvalue_species2\tPvalue_species3\tPvalue_species1&2\tPvalue_species1&3\tPvalue_species2&3\tPvalue_species1&2&3\n");
     for(i=1;i<=NGO;i++){
          if((graph[i-1].exist==1)){
              fprintf(file_out,"GO:%.7lu\t%d",i,graph[i-1].level);
              for (j=1;j<=NSPEC;j++){
                  fprintf(file_out,"\t%lu\t%lu",graph[i-1].count[j-1],graph[i-1].background[j-1]);
              };
              for (j=1;j<=NSTAT;j++){
                  fprintf(file_out,"\t%f",graph[i-1].pvalue[j-1]);
              };
              fprintf(file_out,"\n");
          };    
      };
  } else {
      fprintf(file_out,"Go_term\tLevel\tOccurrences_species1\tBackground_species1\tOccurrences_species2\tBackground_species2\tPvalue_species1\tPvalue_species2\tPvalue_species1&2\n");
      for(i=1;i<=NGO;i++){
          if((graph[i-1].exist==1)){
              fprintf(file_out,"GO:%.7lu\t%d",i,graph[i-1].level);
              for (j=1;j<=2;j++){
                  fprintf(file_out,"\t%lu\t%lu",graph[i-1].count[j-1],graph[i-1].background[j-1]);
              };
              for (j=1;j<=2;j++){
                  fprintf(file_out,"\t%f",graph[i-1].pvalue[j-1]);
              }; j=4; fprintf(file_out,"\t%f",graph[i-1].pvalue[j-1]);
              fprintf(file_out,"\n");
          };    
      };
  };
  
  fclose(file_out);
};

void free_graph(){
  unsigned long j,k;
  
  for(j=1;j<=NGO;j++){
    if (graph[j-1].n_childs>0) free(graph[j-1].childs);
    if (graph[j-1].n_parents>0) free(graph[j-1].parents);
  };    
};
  
//input: go.droso go.danio go.cel bg.go.droso bg.go.danio bg.go.cel output.go.enrichment

int main(int argc, char *argv[])
{
    char *file_set1,*file_set2,*file_set3,*file_bg1,*file_bg2,*file_bg3,*file_output;
    int arg_i, check_i, is_3_species;
    char *init_output="MSGEA_stats.txt";
    
    is_3_species=0;
    check_i=1;
    file_output=init_output;
    for(arg_i=1;arg_i<argc-1;arg_i++)
    {
        if (strcmp(argv[arg_i], "-o") == 0) {arg_i++; file_output=argv[arg_i]; }
        else if (strcmp(argv[arg_i], "-f1") == 0) {arg_i++; check_i=check_i*2; file_set1=argv[arg_i]; }
        else if (strcmp(argv[arg_i], "-f2") == 0) {arg_i++; check_i=check_i*3; file_set2=argv[arg_i]; }
        else if (strcmp(argv[arg_i], "-f3") == 0) {arg_i++; check_i=check_i*11; file_set3=argv[arg_i]; }
        else if (strcmp(argv[arg_i], "-b1") == 0) {arg_i++; check_i=check_i*5; file_bg1=argv[arg_i]; }
        else if (strcmp(argv[arg_i], "-b2") == 0) {arg_i++; check_i=check_i*7; file_bg2=argv[arg_i]; }
        else if (strcmp(argv[arg_i], "-b3") == 0) {arg_i++; check_i=check_i*13; file_bg3=argv[arg_i]; }
    }
    if (check_i==2*3*5*7*11*13) {
        is_3_species=1;
    } else if ((check_i!=2*3*5*7)) {
        fprintf(stderr,"Missing options in command line. Options:\n \t -o output_file : name of output file (default: MSGEA_stats.txt)\n \t -b1 background_file_species1 : file containing background GOs for species 1 (compulsory)\n \t -b2 background_file_species2 : file containing background GOs for species 2 (compulsory)\n \t -b3 background_file_species3 : file containing background GOs for species 3\n \t -f1 file_species1 : file containing selected GOs for species 1 (compulsory)\n \t -f2 file_species2 : file containing selected GOs for species 2 (compulsory)\n \t -f3 file_species3 : file containing selected GOs for species 3\n");
        return(-1);
    }
    printf("Initialize graph of relations between GO terms...\n");
	init_graph();
    printf("Completing initialization...\n");
	init_levels();
	printf("Graph completed.\n");
    printf("Read background GOs for 1st species...\n");
	init_background(file_bg1,1);
    printf("Read background GOs for 2nd species...\n");
	init_background(file_bg2,2);
    if (is_3_species){
        printf("Read background GOs for 3rd species...\n");
        init_background(file_bg3,3);
    };
    printf("Read set GOs for 1st species...\n");
	init_counts(file_set1,1);
    printf("Read set GOs for 2nd species...\n");
    init_counts(file_set2,2);
    if (is_3_species){
        printf("Read set GOs for 3rd species...\n");
        init_counts(file_set3,3);
    };
    printf("Accumulate counts...\n");
	propagate();
    printf("Compute enrichment statistics...\n");
	compute_stat();
    printf("\n");
	printf("Print statistics...\n");
	print_stat(file_output,is_3_species);
	printf("Statistics printed.\n");
	free_graph();
    printf("Run completed.\n");
	return 0;
}
