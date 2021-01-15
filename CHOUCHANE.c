
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

#define Maxx(x, y) (((x) > (y)) ? (x) : (y))
const int INT_MIN = -2147483646;

// ################### TABLO #####################

struct tablo {
  int * tab;
  int size;
};

struct tablo * allocateTablo(int size) {
  struct tablo * tmp = malloc(sizeof(struct tablo));
  tmp->size = size;
  tmp->tab = malloc(size*sizeof(int));
  return tmp;
}

// ################ PRINT ARRAY ##################

void printArray(struct tablo * tmp) {
  int size = tmp->size;
  int i;
  for (i = 0; i < size; ++i) {
    printf("%i ", tmp->tab[i]);
  }
  printf("\n");
}

// ################ REVERSE ARRAY ##################
 // https://www.programmingsimplified.com/c-program-reverse-array
void reverse_array(int a[], int n)
{
  int c, t;

  for (c = n/2; c < n/4; c++) {
    t = a[c];                  
    a[c] = a[n-c-1];
    a[n-c-1] = t;
  }
}

// ################### MONTEE #####################

void montee(struct tablo * source, struct tablo * destination) {

  int sourceSize = source->size;
  int m = log10(sourceSize)/log10(2);

  #pragma omp parallel for
  for(int i= 0; i<sourceSize;i++){
    destination->tab[sourceSize+i] = source->tab[i];
  }
  
  //FIXME provoque des warning Valgrind
  // printArray(destination);

  for (int i=m - 1 ; i >= 0 ; i--) {

    int Max1 = pow(2,i+1)-1;

    #pragma omp parallel for
    for (int j=pow(2,i) ;j <= Max1 ; j++) {
      destination->tab[j] = destination->tab[2*j] + destination->tab[2*j+1];	
    }
  }       
}

void monteeMax(struct tablo * source, struct tablo * destination) {

  int sourceSize = source->size;
  int m = log10(sourceSize)/log10(2);

  #pragma omp parallel for
  for(int i= 0; i<sourceSize;i++){
    destination->tab[sourceSize+i] = source->tab[i];
  }
  
  //FIXME provoque des warning Valgrind
  // printArray(destination);

  for (int i=m - 1 ; i >= 0 ; i--) {

    int Max1 = pow(2,i+1)-1;

    #pragma omp parallel for
    for (int j=pow(2,i) ;j <= Max1 ; j++) {
      destination->tab[j] = Maxx(destination->tab[2*j] , destination->tab[2*j+1]);	
    }
  }       
}

// ################### DESCENTE #####################

void descente(struct tablo * a, struct tablo * b) {

  int size = a->size/2;

  //FIXME : manque l élement neutre
  int m = log10(size)/log10(2);
  b->tab[1]=0;

  for (int i= 1; i <= m ; i++) {
    
    int Max2 = (pow(2,i+1)-1);
    #pragma omp parallel for
    for (int j= pow(2,i) ; j <= Max2 ; j++) {
      
      if( j%2 == 0) {
        b->tab[j] = b->tab[j/2];
      } 
      else {
        b->tab[j] = b->tab[(j-1)/2] + a->tab[j-1];
      }
    }
  }
}

void descenteMax(struct tablo * a, struct tablo * b) {

  int size = a->size/2;

  //FIXME : manque l élement neutre
  int m = log10(size)/log10(2);
  b->tab[1]=INT_MIN;

  for (int i= 1; i <= m ; i++) {
    
    int Max2 = (pow(2,i+1)-1);
    #pragma omp parallel for
    for (int j= pow(2,i) ; j <= Max2 ; j++) {
      
      if( j%2 == 0) {
        b->tab[j] = b->tab[j/2];
      } 
      else {
        b->tab[j] = Maxx(b->tab[(j-1)/2] , a->tab[j-1]);
      }
    }
  }
}

// ################### FINAL #####################

void final(struct tablo * a, struct tablo *b) {
  
  int size = a->size;
  int m = log10(size)/log10(2);
  
  int Max3 = (pow(2,m)-1);
  #pragma omp parallel for
  for(int j = pow(2,m-1) ; j <= Max3 ; j++ )
  {
    b->tab[j] = b->tab[j] + a->tab[j];
  }
}

void finalMax(struct tablo * a, struct tablo *b) {
  
  int size = a->size;
  int m = log10(size)/log10(2);
  
  int Max3 = (pow(2,m)-1);
  #pragma omp parallel for
  for(int j = pow(2,m-1) ; j <= Max3 ; j++ )
  {
    b->tab[j] = Maxx(b->tab[j] , a->tab[j]);
  }
}

// ################### GENERATE TEST #####################

void generateArray(struct tablo * s, struct tablo * c) {
  //construction d'un tableau pour tester
  s->size=16;
  s->tab=malloc(s->size*sizeof(int));
  c->size=16;
  c->tab=malloc(s->size*sizeof(int));

  s->tab[0]=3;  c->tab[15]=3;
  s->tab[1]=2;  c->tab[14]=2;
  s->tab[2]=-7; c->tab[13]=-7;
  s->tab[3]=11; c->tab[12]=11;
  s->tab[4]=10; c->tab[11]=10;
  s->tab[5]=-6; c->tab[10]=-6;
  s->tab[6]=4;  c->tab[9]=4;
  s->tab[7]=9;  c->tab[8]=9;
  s->tab[8]=-6; c->tab[7]=-6;
  s->tab[9]=1;  c->tab[6]=1;
  s->tab[10]=-2;c->tab[5]=-2;
  s->tab[11]=-3;c->tab[4]=-3;
  s->tab[12]=4; c->tab[3]=4;
  s->tab[13]=-3;c->tab[2]=-3;
  s->tab[14]=0; c->tab[1]=0;
  s->tab[15]=2; c->tab[0]=2;
}

// ################### PREFIX SUM #####################

void prefixSum(struct tablo * in, struct tablo *out){

  struct tablo * tmp = malloc(sizeof(struct tablo));
  tmp->tab = malloc(in->size*2*sizeof(int));
  tmp->size =in->size*2;
  tmp->tab[0]=0;

  struct tablo * preout = malloc(sizeof(struct tablo));
  preout->tab = malloc(in->size*2*sizeof(int));
  preout->size =in->size*2;
  preout->tab[0]=0;

  montee(in , tmp);

  descente(tmp, preout);	

  final(tmp,preout);
  
  #pragma omp parallel for
  for (int i = 0; i < in->size; ++i)
  {
  	out->tab[i] = preout->tab[in->size+i];
  }

  free(tmp->tab);
  free(tmp);

}

// ################### SUFIX SUM #####################

void sufixSum(struct tablo * in, struct tablo *out){

  struct tablo * tmp = malloc(sizeof(struct tablo));
  tmp->tab = malloc(in->size*2*sizeof(int));
  tmp->size =in->size*2;
  tmp->tab[0]=0;

  struct tablo * preout = malloc(sizeof(struct tablo));
  preout->tab = malloc(in->size*2*sizeof(int));
  preout->size =in->size*2;
  preout->tab[0]=0;

  montee(in , tmp); //***

  descente(tmp, preout);

  final(tmp,preout);

  #pragma omp parallel for
  for (int i = 0; i < in->size; ++i)
  {
  	out->tab[i] = preout->tab[in->size*2-1-i];
  }

  free(tmp->tab);
  free(tmp);
  free(preout->tab);
  free(preout);
}

// ################### PREFIX MAX #####################

void prefixMax(struct tablo * in, struct tablo *out){

  struct tablo * tmp = malloc(sizeof(struct tablo));
  tmp->tab = malloc(in->size*2*sizeof(int));
  tmp->size =in->size*2;
  tmp->tab[0]=0;

  struct tablo * preout = malloc(sizeof(struct tablo));
  preout->tab = malloc(in->size*2*sizeof(int));
  preout->size =in->size*2;
  preout->tab[0]=0;

  monteeMax(in , tmp);

  descenteMax(tmp, preout);	

  finalMax(tmp,preout);

  #pragma omp parallel for
  for (int i = 0; i < in->size; ++i)
  {
  	out->tab[i] = preout->tab[in->size+i];
  }

  free(tmp->tab);
  free(tmp);

}

// ################### SUFIX MAX #####################

void sufixMax(struct tablo * in, struct tablo *out){

  struct tablo * tmp = malloc(sizeof(struct tablo));
  tmp->tab = malloc(in->size*2*sizeof(int));
  tmp->size =in->size*2;
  tmp->tab[0]=0;

  struct tablo * preout = malloc(sizeof(struct tablo));
  preout->tab = malloc(in->size*2*sizeof(int));
  preout->size =in->size*2;
  preout->tab[0]=0;

  monteeMax(in , tmp); //***

  descenteMax(tmp, preout);

  finalMax(tmp,preout);

  #pragma omp parallel for
  for (int i = 0; i < in->size; ++i)
  {
  	out->tab[i] = preout->tab[in->size*2-1-i];
  }

  free(tmp->tab);
  free(tmp);
  free(preout->tab);
  free(preout);
}

// ############### LARGEST SUBARRAY ##################

//working on this part still not functional

	void maxSubArraySumPrint(int a[],struct tablo * Q, int size) 
	{ 
   int maxM_Value= 0, maxValue = 0;
   int nb_maxValue = 0;

   for (int i = 0; i < size; ++i)
   {
   	if (a[i]>maxM_Value) 
   		{
   			maxM_Value = a[i];
   			nb_maxValue = 1; 
   		}
   		else if (a[i]== maxM_Value)  
   			nb_maxValue = nb_maxValue+1;
   }

   int mAxValArray[nb_maxValue];

   for (int i = 0; i < size; ++i)
   {   	
   	if (a[i] == maxM_Value)
   	{
   		maxValue = maxValue + Q->tab[a[i]];
   		mAxValArray[i] =  Q->tab[a[i]];
   	}
   }


	printf("%i",maxValue);
	for (int i = 0; i < nb_maxValue; ++i)
	{
		printf(" %i",mAxValArray[i]);
	}
	printf("\n");
	} 

// ##################### MAIN ########################

int main(int argc, char **argv) {
  struct tablo source, reverseSource;
  generateArray(&source,&reverseSource);
 
  int MS[source.size], MP[source.size], M[source.size];



  struct tablo * PSUM = malloc(sizeof(struct tablo));
  PSUM->tab = malloc(source.size*sizeof(int));
  PSUM->size =source.size;
  PSUM->tab[0]=0;

  struct tablo * SSUM = malloc(sizeof(struct tablo));
  SSUM->tab= malloc(source.size*sizeof(int));
  SSUM->size=source.size;
  SSUM->tab[0]=0;

    struct tablo * SMAX = malloc(sizeof(struct tablo));
  SMAX->tab = malloc(source.size*sizeof(int));
  SMAX->size =source.size;
  SMAX->tab[0]=0;

  struct tablo * PMAX = malloc(sizeof(struct tablo));
  PMAX->tab= malloc(source.size*sizeof(int));
  PMAX->size=source.size;
  PMAX->tab[0]=0;

  prefixSum(&source, PSUM);

  sufixSum(&reverseSource, SSUM);

  sufixMax(PSUM, SMAX);

  prefixMax(SSUM, PMAX);


   #pragma omp parallel for
   for (int i = 0; i < source.size; ++i)
   {
   	MS[i] = PMAX->tab[i] - SSUM->tab[i] + source.tab[i];
   	MP[i] = SMAX->tab[i] - PSUM->tab[i] + source.tab[i];
   	M[i]  = MS[i] + MP[i] - source.tab[i];
   }

   maxSubArraySumPrint(M,&source, source.size);


free(PSUM->tab);
free(PSUM);
free(SSUM->tab);
free(SSUM);
free(SMAX->tab);
free(SMAX);
free(PMAX->tab);
free(PMAX);
}
