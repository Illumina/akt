#include <assert.h>
#include "pedigree.h"

//#define DEBUG_PEDIGREE

sampleInfo::sampleInfo(string fname)
{
    readFromPlinkFam(fname);
    buildIndex();
}


sampleInfo::sampleInfo(string fname,bcf_hdr_t *hdr)
{
    readFromPlinkFam(fname);

    string sample_list=id[0];
    for(vector<string>::iterator it=id.begin()+1;it!=id.end();it++)
    {
        sample_list += ","+(*it);
    }

    int new_N=bcf_hdr_nsamples(hdr);
    vector<string> new_fid(new_N,""),new_id(new_N,""),new_dad(new_N,""),new_mum(new_N,"");
    vector<int> new_sex(new_N,3),new_status(new_N,-9);
    for(int i=0;i<new_N;i++)
    {
        new_id[i]=(string)hdr->samples[i];
        int idx=0;
        while(idx<N && id[idx]!=new_id[i])
        {
            idx++;
        }
        if(idx<N)
        {
            new_fid[i] = fid[idx];
            new_dad[i] = dad[idx];
            new_mum[i] = mum[idx];
            new_sex[i] = sex[idx];
            new_status[i] = status[idx];
        }
        //else sample wasnt in original pedigree
    }
    sex=new_sex;
    status=new_status;
    id=new_id;
    fid=new_fid;
    dad=new_dad;
    mum=new_mum;
    N=new_N;
    cerr << "Found "<<N<<" in both the pedigree and bcf."<<endl;
    buildIndex();
}

int sampleInfo::readFromPlinkFam(string fname)
{
    std::ifstream inf(fname.c_str(), std::ifstream::in);
    if(!inf)
    {
        cerr << "ERROR: problem reading "<<fname<<endl;
        exit(1);
    }
    N = 0;
    string line;
    stringstream ss;
    vector<string> entry;
    while(getline(inf,line))
    {
//        cerr << N << " " <<line <<endl;
        stringSplit(line,entry);
//        for(size_t  i=0;i<entry.size();i++) cerr << i << " " <<entry[i]<<endl;
        assert(entry.size()>=4);
        fid.push_back(entry[0]);
        id.push_back(entry[1]);
        dad.push_back(entry[2]);
        mum.push_back(entry[3]);
        sex.push_back(atoi(entry[4].c_str()));
        if(entry.size()>=6)
        {
            status.push_back(atoi(entry[5].c_str()));
        }
        else
        {
            status.push_back(-9);
        }
        N++;
    }
    cerr << "Read "<<N<<" individuals from "<<fname<<endl;
    if(N==0) die("no samples read from "+fname);
    return (0);
}

int sampleInfo::buildIndex()
{
    dadidx.assign(N,-1);
    mumidx.assign(N,-1);
    nduo=ntrio=0;
    for(int i=0;i<N;i++)
    {
        if( dad[i]!="0" )
        {
            dadidx[i] = 0;
            while(dadidx[i]<N&&dad[i]!=id[dadidx[i]])
            {
                dadidx[i]++;
            }
            if(dadidx[i]>=N)
            {
                dadidx[i]=-1;
            }
        }

        if( mum[i]!="0" )
        {
            mumidx[i] = 0;
            while(mumidx[i]<N&&mum[i]!=id[mumidx[i]])
            {
                mumidx[i]++;
            }
            if(mumidx[i]>=N)
            {
                mumidx[i]=-1;
            }
        }

        if(mumidx[i]!=-1&&dadidx[i]!=-1)
        {
            ntrio++;
        }
        else if(mumidx[i]!=-1||dadidx[i]!=-1)
        {
            nduo++;
        }
    }
    cerr<<"Found "<<ntrio<<" trios and "<<nduo<<" duos"<<endl;

    for(int i=0;i<N;i++)
    {
        if(mumidx[i]!=-1 && dadidx[i]!=-1)
        {
            pair<int,int>key(dadidx[i],mumidx[i]);
            if(!parent_map.count(key))
            {
                parent_map[key] = vector<int>();
            }
            parent_map[key].push_back(i);
        }
    }
    int num_affected=0;
    int num_unaffected=0;
    int num_nuclear=0;
    int num_sample_in_nuc=0;
    for( map<pair<int,int>, vector<int>  >::iterator it1=parent_map.begin();it1!=parent_map.end();it1++)
    {
#ifdef DEBUG_PEDIGREE
        cerr << id[it1->first.first] << " " << id[it1->first.second] << " |";
      for(size_t i=0;i<it1->second.size();i++)
      {
	  cerr << " " <<id[it1->second[i]];
      }	  
      cerr<<endl;
#endif
        for(size_t i=0;i<it1->second.size();i++)
        {
            if(status[it1->second[i]]==0) num_affected++;
            if(status[it1->second[i]]==2) num_unaffected++;
        }
        num_nuclear++;
        num_sample_in_nuc+=2+it1->second.size();
    }

    cerr << num_nuclear << " nuclear families (both parents assayed) containing "<<num_sample_in_nuc<<" individuals"<<endl;
    cerr << num_affected << " affected children"<<endl;
    cerr << num_unaffected << " unaffected children"<<endl;
    return(0);
}

int sampleInfo::getMumIndex(int idx)
{
    if (mumidx[idx] != -1)
    {
        return(mumidx[idx]);
    }
    else
    {
        return(-1);
    }
}

int sampleInfo::getDadIndex(int idx)
{
    if (dadidx[idx] != -1)
    {
        return(dadidx[idx]);
    }
    else
    {
        return(-1);
    }
}

