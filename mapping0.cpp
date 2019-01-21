#include <iostream>
#include <fstream>
#include <algorithm>//std::max_element
#include <string> //std::max_element
#include <cstring>
#include <stdio.h>
#include<regex>
#include <set>
//#include <limits>

using namespace std;

int graph_alignment( char *Protein, string kmer_array[], int node_array[][4], int I, int J, int k){
  int M[I+1][J+1];
  int delta_match = 2;
  int delta_mismatch = -1;
  int delta = 0;
  M[0][0] = 0;
  for (int i=0; i<I+1; i++){
      M[i][0] = 0;
  }
  for (int j=0; j<J+1; j++){
      M[0][j] = 0;
  }

  for (int i=1; i<I+1; i++){
    for (int nodej=1; nodej < J+1; nodej++){
      if (Protein[i-1] == kmer_array[nodej-1][k-1]){
        delta = delta_match;
      }
      else{
        delta = delta_mismatch;
      }
      //cout << " " << kmer_array[nodej-1] << endl;
      int pre_num = 1;
      for (int prei=0;prei<4;prei++){
        if (node_array[nodej-1][prei] > -1){
          pre_num++;
        }
      }
      
      int myints[pre_num+1];
      int p=0;
      for (int prei=0;prei<4;prei++){
        if (node_array[nodej-1][prei] > -1){
          myints[p] = M[i-1][node_array[nodej-1][prei]+1] + delta;
          p++;
        }
      }
      myints[p] = 0;
      M[i][nodej] = *max_element(myints, myints + pre_num);
    }
  }

  cout << "dynamic programming table: " << endl;
  int largest_score = 0;
  int largest_i = I, largest_j = J;
  cout <<"m matrix: " << endl;
    for (int i=0; i<I+1; i++){
      for (int nodej=0; nodej<J+1; nodej++){
        if (M[i][nodej] >= largest_score){
          largest_score = M[i][nodej];
          largest_i = i;
          largest_j = nodej;
        }
        cout << M[i][nodej];
        cout <<"\t| ";
      }
      cout <<endl;
    }
  cout <<"end table" << endl;
  cout <<endl;

  cout << "alignment sequence:" << endl;
  vector<char> protein_list;
  vector<char> dna_list;
  protein_list.push_back(' ');
  dna_list.push_back(' ');
  int i = largest_i, nodej = largest_j;
  int track_score = largest_score;
  while (M[i][nodej] != 0){
    cout << "M[" << i << "][" << nodej << "] = "<< M[i][nodej] << endl;
    for (int pre_j=0;pre_j<4;pre_j++){
      if (node_array[nodej-1][pre_j] > -1){
        if (track_score == M[i][nodej]){
          if ((M[i][nodej] == M[i-1][node_array[nodej-1][pre_j]+1] + delta_match) && (Protein[i-1] == kmer_array[nodej-1][k-1])){
            protein_list.push_back(Protein[i-1]);
            dna_list.push_back(kmer_array[nodej-1][k-1]);
            i--;
            nodej = node_array[nodej-1][pre_j]+1;
            track_score = M[i][nodej];
            break;
          }
          else if ((M[i][nodej] == M[i-1][node_array[nodej-1][pre_j]+1] + delta_mismatch) && (Protein[i-1] != kmer_array[nodej-1][k-1])){
            protein_list.push_back(Protein[i-1]);
            dna_list.push_back(kmer_array[nodej-1][k-1]);
            i--;
            nodej = node_array[nodej-1][pre_j]+1;
            track_score = M[i][nodej];
            break;
          }
        }
      }
    }
  }
  dna_list.push_back(' ');
  protein_list.push_back(' ');

  cout << endl;
  for (std::vector<char>::const_iterator i = dna_list.end()-1; i != dna_list.begin(); i--){
    cout << *i ;
  }
  cout << endl;
  for (std::vector<char>::const_iterator i = protein_list.end()-1; i != protein_list.begin(); i--){
    cout << *i ;
  }
  cout << endl;

  return largest_score;
}

int binarySearch(string names[], int size, string value){
  int first = 0,
      last = size - 1,
      middle,
      position = -1;
  bool found = false;

  while (!found && first <= last){
    middle = (first + last) / 2;
    if (names[middle] == value){
        found = true;
        position = middle;
    }
    else if (names[middle] > value){
        last = middle - 1;
    }
    else {
        first = middle + 1;
    }
  }
  return position;
}

int main (int argc,char *argv[]) {
  if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <DNA sequence> <Kmer File>" << std::endl;
        return 1;
    }

  ifstream dnasequence (argv[1]);
  ifstream kmerfile (argv[2]);
  if (!dnasequence.is_open())
    cout<<"Could not open protein file\n";
  if (!kmerfile.is_open())
    cout<<"Could not open graph file\n";

  std::string line1;
  getline(dnasequence,line1); 
  int m = line1.length(); 
  char Protein[m+1]; 
  strcpy(Protein, line1.c_str());

  string kmer;
  int kmer_count = 0;
  int kmer_length;

  while (kmerfile >> kmer){
    if (regex_match(kmer, regex("(>)(.*)"))){
      //cout << i << endl;
    }
    else{
      kmer_count++;
      kmer_length = kmer.length();
    }
  }
  string kmer_list[kmer_count];
  string key_list[kmer_count];
  kmerfile.clear();
  kmerfile.seekg(0, kmerfile.beg);

  int i=0;
  while (kmerfile >> kmer){
    if (regex_match(kmer, regex("(>)(.*)"))){
      //cout << i << endl;
    }
    else{
      kmer_list[i]=kmer;
      string kmer_postk = kmer.substr(1,kmer_length-1);
      string firstletter = kmer.substr(0,1);
      string key = kmer_postk.append(firstletter);
      key_list[i]=key;
      i++;
    }
  }
  sort(key_list, key_list + kmer_count);
  int node_array[kmer_count][4];
  string kmer_array[kmer_count];

  //binary search
  for (int j=0;j<kmer_count;j++){
    string tmpnode = kmer_list[j].substr(1,kmer_length-1);
    tmpnode.append(kmer_list[j].substr(0,1));
    int node_i = binarySearch(key_list, kmer_count, tmpnode);
    kmer_array[node_i] = kmer_list[j];
  }
  for (int j=0;j<kmer_count;j++){
    string tmpstrA = kmer_array[j].substr(0,kmer_length-1).append("A");
    node_array[j][0] = binarySearch(key_list, kmer_count, tmpstrA);

    string tmpstrC = kmer_array[j].substr(0,kmer_length-1).append("C");
    node_array[j][1] = binarySearch(key_list, kmer_count, tmpstrC);

    string tmpstrG = kmer_array[j].substr(0,kmer_length-1).append("G");
    node_array[j][2] = binarySearch(key_list, kmer_count, tmpstrG);

    string tmpstrT = kmer_array[j].substr(0,kmer_length-1).append("T");
    node_array[j][3] = binarySearch(key_list, kmer_count, tmpstrT);
  }

  //cout << "Graph Node: " << endl;
  //for (int i=0;i<n;i++){
  //  cout << "PreNodes: " << endl;
  //  for (int j=0;j<4;j++){
  //    cout << node_array[i][j] << " ";
  //  }
  //  cout << endl;
  //}
  cout << "Global alignment score is " << graph_alignment( Protein, kmer_array, node_array, m, kmer_count, kmer_length) << endl;
  cout << endl;
  return 0;
}
