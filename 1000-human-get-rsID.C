/*
input:
(1) a file name (which lists diff chrom VCF files)
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099
		16	53842908	rs12149832	G	A	100	PASS	AC=1143;AF=0.228235;AN=5008;NS=2504;DP=20175;EAS_AF=0.1706;AMR_AF=0.232;AFR_AF=0.0582;EUR_AF=0.4185;SAS_AF=0.319;AA=G|||;VT=SNP	GT	1|0	1|1	1|0

(2) individual ID - poupulation file (HG00096	europe)
(3) obes rs ID file (rs1011731	G)
(4) # of fields in header/data 

output file:
rs1011731	G(GWAS allele)	A(alt allele)	QUAL    FILTER  AFR	100(G allele count)	300(A count)	400(total)	1/4(G freq)	AMR	200	200	400	1/2	ASN	..	EUR	..	SAN	..	ALL	2504(G)	2504(A)	50008	1/2

rsXXXX not found in 1000-genome genotype table
rsXXXX G allele not found in 1000-genome genotype table


map<string, char> snpMap //<"rs1011731"  'G'>
map<string, bool> snpFound //<"rs1011731"  FALSE>

read obes rs ID file, store snp in the snpMap

map<string, string> idPopu <"HG00096", "AFR">
set<string> popuSet
struct twoInt {int gwasAllelFreq; int altAllelFreq};
map<string, twoInt> allelFreq; //<"AFR", {100, 300}>

read individual ID - poupulation file
  populate idPopu
  populate popuSet



vector<string> vcfHeader;
vector<string> vcfFiles   // store all vcf files
map<string, twoInt> allelFreq; //<"AFR", {100, 300}>

header flag 'false'
open a vcf file
 read a line
  if header flag 'false'
    if the line is a header
     header flag 'true'
     vcfHeader-push_back(); a[20]="HG00096"
  else //already encounted the header line
   if  the line contains rs ID, if the snp allele is in the genotype, 
      snpFound set TRUE
      for each population, twoInt cleared to {0,0}
      a[20]="0|1" a[3]="ref allele G:0"  a[4]="alt allele A:1"
      vcfHeader[20]="HG00096", idPopu["HG00096"]="AFR", b='0' c='1'
      if(snp allele same as ref allele)
        if(b == '0')
          allelFreq["AFR"].gwasAllelFreq++;
        if(b == '1')
          allelFreq["AFR"].altAllelFreq++;
        if(c ==
      output population allele results  

close a vcf file

output rs ID not found in the VCF files

 */

#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <set>
#include <vector>
#include <sstream>
#include <map>

using namespace std;

struct twoInt {int gwasAllelFreq;
                int altAllelFreq;};


void readFile(const char * argv, ifstream & input);
void writeFile(const char * argv, ofstream & output);
void getFieldContent(vector<string> & fields, char del, string line);
void getFieldContent2(vector<string> & fields, string del, string line);

int main(int argc, char *argv[])
{
 if(argc != 4)
 {
  cerr << argv[0] << "   1[a file name (which lists necessary diff chrom VCF files)]   2[obes rs ID file (rs1011731)]  3[output-file]" << endl;
  cerr << "******************************************************************************" << endl;
  cerr<<"warning: all VCF files are assumed to have same format" << endl;
  cerr << "Format (order and number of fields are assumed): " << endl;
  cerr<<"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099	HG00100	HG00101	..." << endl;
  return 1;
 }

 ifstream input;
 ofstream output;
 string line;
 vector<string> lineFields;
 map<string, char> snpMap; //<"rs1011731"  'G'>
 map<string, bool> snpFound; //<"rs1011731"  FALSE>
 map<string, string> idPopu; // <"HG00096", "AFR">
 set<string> popuSet;
 map<string, twoInt> allelFreq; //<"AFR", {100, 300}>
 vector<string> vcfHeader;
 vector<string> vcfFiles;   // store all vcf files
 bool vcfHeaderBool; 

 readFile(argv[2], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent2(lineFields, " \n\t", line);
  if(lineFields.size() == 1)
  {
   snpMap[lineFields[0] ] = '@';
   snpFound[lineFields[0] ] = false;
  }
  getline(input, line);
 }
 input.close();
 cout << snpMap.size() << " unique rsIDs." << endl;


 readFile(argv[1], input);
 getline(input, line);
 while(!input.eof())
 {
  getFieldContent2(lineFields, "\n\t", line);
  for(int i = 0; i < lineFields.size(); i++)
  {
   vcfFiles.push_back(lineFields[i]);
  }
  getline(input, line);
 }
 input.close();
 cout << vcfFiles.size() << " chromosome vcf files." << endl;

 writeFile(argv[3], output);
 vcfHeaderBool = false;
 for(int i = 0; i < vcfFiles.size(); i++)
 {
  readFile(vcfFiles[i].c_str(), input);
  getline(input, line);
  while(!input.eof())
  {
   getFieldContent(lineFields, '\t', line);
   if(lineFields.size() >= 3 &&  snpMap.count(lineFields[2]) > 0 )
   {
    snpFound[lineFields[2]] = true;
    output << line << endl;
   }
   getline(input, line);
  }
  input.close();
 }
 for(map<string, bool> ::iterator it = snpFound.begin(); it != snpFound.end(); it++)
 {
  if(it->second == false)
  {
   cerr << it->first << " not found in 1000-genome genotype table" << endl;
  }
 }
 output.close();
 return 0;
}

void readFile(const char * argv, ifstream & input)
{
 input.open(argv, ios::in);
 if(!input)
 {
  cerr << argv << " can not be opened for reading.\n";
  exit(1);
 }
}

void writeFile(const char * argv, ofstream & output)
{
 output.open(argv, ios::out);
 if(!output)
 {
  cerr << argv << " can not be opened for writing.\n";
  exit(1);
 }
}

void getFieldContent(vector<string> & fields, char del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del == line[i])
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   str = line.substr(pos[i] + 1, len);
   fields.push_back(str);
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}

void getFieldContent2(vector<string> & fields, string del, string line)
{
 vector<int> pos;
 string str;
 int len, size;

 fields.clear();
 for(int i = 0; i < line.length(); i++)
 {
  if(del.find(line[i]) != string::npos)
    pos.push_back(i);
 }
 if(pos.size() > 0)
 {
  len = pos[0];
  str = line.substr(0, len);
  if(len > 0)
    fields.push_back(str);
  for(int i = 0; i < pos.size() - 1; i++)
  {
   len = pos[i+1] - pos[i] - 1;
   if(len > 0)
   {
    str = line.substr(pos[i] + 1, len);
    fields.push_back(str);
   }
  }
  size = pos.size();
  if(pos[size-1] < line.length() - 1) //not at the end of line
  {
   str = line.substr(pos[size-1] + 1);
   fields.push_back(str);
  }
 }
 else
 {
  fields.push_back(line);
 }
}


