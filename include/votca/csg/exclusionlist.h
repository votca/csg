/* 
 * Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _exclusionlist_H
#define	_exclusionlist_H

#include <memory>
#include <iostream>
#include <list>
#include <map>
#include "bead.h"

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

/// \todo fill _excl_by_bead
/// \todo no ids but pointers, use PairList

class Topology;
class Bead;

class ExclusionList
{
public:
    ExclusionList() {}
    ~ExclusionList() { Clear(); }
    
    void Clear(void);

    template<typename iteratable>
    void Remove(iteratable &l);

    template<typename iteratable>
	void ExcludeList(iteratable &l);
    
    struct exclusion_t {
        std::shared_ptr<Bead> _atom;
        list<std::shared_ptr<Bead>> _exclude;
    };    

    void CreateExclusions(Topology *top);
    exclusion_t *GetExclusions(std::shared_ptr<Bead> bead);
    
    typedef list< exclusion_t * >::iterator iterator;
    
    iterator begin() { return _exclusions.begin(); }
    iterator end() { return _exclusions.end(); }
    
    bool IsExcluded(std::shared_ptr<Bead> bead1,std::shared_ptr<Bead> bead2);

    template<typename iteratable>
    void InsertExclusion(std::shared_ptr<Bead> bead, iteratable &excluded);

    void InsertExclusion(std::shared_ptr<Bead> bead1, std::shared_ptr<Bead> bead2);

    void RemoveExclusion(std::shared_ptr<Bead> bead1, std::shared_ptr<Bead> bead2);
private:
    list< exclusion_t * > _exclusions;
    map<std::shared_ptr<Bead> , exclusion_t *> _excl_by_bead;
    
    friend std::ostream &operator<<(std::ostream &out, ExclusionList& exl);
};

inline ExclusionList::exclusion_t * ExclusionList::GetExclusions(std::shared_ptr<Bead> bead)
{
   map<std::shared_ptr<Bead> , exclusion_t *>::iterator iter  = _excl_by_bead.find(bead);
   if(iter == _excl_by_bead.end()) return NULL;
   return (*iter).second;
}

template<typename iteratable>
inline void ExclusionList::Remove(iteratable &l)
{
    typename iteratable::iterator i, j;

    for ( i = l.begin(); i != l.end(); ++i ) {
    	for ( j = i; j != l.end(); ++j ) {
    		RemoveExclusion(*i, *j);
    	}
    }
}

template<typename iteratable>
inline void ExclusionList::ExcludeList( iteratable &l ) {
    typename iteratable::iterator i, j;

    for ( i = l.begin(); i != l.end(); ++i ) {
    	for ( j = i; j != l.end(); ++j ) {
    		InsertExclusion(*i, *j);
    	}
    }
}

template<typename iteratable>
inline void ExclusionList::InsertExclusion(std::shared_ptr<Bead> bead1_, iteratable &l)
{
	for(typename iteratable::iterator i=l.begin(); i!=l.end(); ++i) {
    std::shared_ptr<Bead> bead1 = bead1_;
    std::shared_ptr<Bead> bead2 = *i;
		if (bead2->getId() < bead1->getId()) swap(bead1, bead2);
		if(bead1==bead2) continue;
		if(IsExcluded(bead1, bead2)) continue;
		exclusion_t *e;
		if((e = GetExclusions(bead1)) == NULL) {
			e = new exclusion_t;
			e->_atom = bead1;
			_exclusions.push_back(e);
			_excl_by_bead[ bead1 ] = e;
		}
		e->_exclude.push_back(bead2);
	}
}

//template<>
inline void ExclusionList::InsertExclusion(std::shared_ptr<Bead> bead1,std::shared_ptr<Bead> bead2) {
    if (bead2->getId() < bead1->getId()) swap(bead1, bead2);
	if(bead1==bead2) return;
	if(IsExcluded(bead1, bead2)) return;

	exclusion_t *e;
	if((e = GetExclusions(bead1)) == NULL) {
		e = new exclusion_t;
		e->_atom = bead1;
		_exclusions.push_back(e);
		_excl_by_bead[ bead1 ] = e;
	}
	e->_exclude.push_back(bead2);
}

inline void ExclusionList::RemoveExclusion(std::shared_ptr<Bead> bead1,std::shared_ptr<Bead> bead2) {
    if (bead2->getId() < bead1->getId()) swap(bead1, bead2);
    if(bead1==bead2) return;
    if(!IsExcluded(bead1, bead2)) return;
    list<exclusion_t*>::iterator ex;
    for(ex=_exclusions.begin(); ex!=_exclusions.end(); ++ex)
        if((*ex)->_atom == bead1) break;
    if(ex==_exclusions.end()) return;
    (*ex)->_exclude.remove(bead2);
    if((*ex)->_exclude.empty()) {
        (*ex)=NULL;
        _exclusions.erase(ex);
    }
    _exclusions.remove(NULL);
}

std::ostream &operator<<(std::ostream &out,ExclusionList& ex);

}}

#endif	/* _exclusionlist_H */



