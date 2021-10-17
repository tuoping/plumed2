/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2021 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "VesselRegister.h"
#include "Vessel.h"
#include "StoreDataVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

class Contents : public Vessel {
private:
  // unsigned mycomponent;
  StoreDataVessel* mystash;
  std::vector<unsigned> members;
  std::vector<Value*> value_out;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit Contents( const vesselbase::VesselOptions& da );
  std::string description() override {return "Contents";};
  void resize() override {};
  void calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const override {}
  void finish( const std::vector<double>& buffer ) override;
  bool applyForce( std::vector<double>& forces ) override {return false; };
};

PLUMED_REGISTER_VESSEL(Contents,"CONTENTS")

void Contents::registerKeywords( Keywords& keys ) {
  Vessel::registerKeywords( keys );
  keys.remove("LABEL");
  keys.add("compulsory","CONTENTS","the list of contents that you would like to calculate");
}

void Contents::reserveKeyword( Keywords& keys ) {
  keys.reserve("optional","CONTENTS","calculate the contents of collective variables. "
               "If you would like the second and third contents of the third component you would use CONTENTS={CONTENTS=2-3}.  The contents would then be referred to "
               "using the labels content-3 and content-3.  This syntax is also required if you are using numbered CONTENT keywords i.e. CONTENTS1, CONTENTS2...");
  keys.reset_style("CONTENTS","vessel");
  keys.addOutputComponent("content","CONTENTS","the central contents of the distribution of values. The second content "
                          "would be referenced elsewhere in the input file using "
                          "<em>label</em>.content-2, the third as <em>label</em>.content-3, etc.");
}


Contents::Contents( const vesselbase::VesselOptions& da) :
  Vessel(da)
{
  mystash = getAction()->buildDataStashes( NULL );
  ActionWithValue* a=dynamic_cast<ActionWithValue*>( getAction() );
  plumed_massert(a,"cannot create passable values as base action does not inherit from ActionWithValue");

  std::vector<std::string> contents; std::string valstr;

  valstr = "content-";
  if( getNumericalLabel()==0 ) {
    contents = Tools::getWords(getAllInput(), "\t\n ,");
    Tools::interpretRanges(contents);// mycomponent=1;
  } else {
    parseVector("CONTENTS", contents);
    Tools::interpretRanges(contents);
  }
  unsigned nn;
  for(unsigned m=0; m<contents.size(); ++m) {
    for(unsigned i=0; i<mystash->getNumberOfStoredValues(); ++i) {
      a->addComponent( valstr + contents[i] );
      value_out.push_back( a->copyOutput( a->getNumberOfComponents()-1 ) );    
    }
    Tools::convert( contents[m], nn );
    members.push_back( nn ); std::string num; Tools::convert(members[m],num);
  }
}


void Contents::finish( const std::vector<double>& buffer ) {
  std::vector<double>  myvalues( getAction()->getNumberOfQuantities() );

  Value myvalue;
  if(getAction()->isPeriodic()){
      std::string str_min, str_max;
      getAction()->retrieveDomain(str_min, str_max);
      double min, max;
      Tools::convert(str_min, min);
      Tools::convert(str_max, max);
      myvalue.setDomain(str_min, str_max);
  } else {
      myvalue.setNotPeriodic();
  }

  unsigned k=0;
  unsigned nn;
  for(unsigned m=0; m<members.size(); ++m) {
    nn = members[m];
    for(unsigned i=0; i<mystash->getNumberOfStoredValues(); ++i) {
      mystash->retrieveSequentialValue( i, false, myvalues );
      double tmp=myvalues[nn];
      value_out[k]->set(tmp);
      k++;
    }
  }
}


}
}
