#ifndef __PYTHIAEVENT_H__
#define __PYTHIAEVENT_H__

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"



//----------------------------------------------------------------------------------------
class PythiaParticle : public TObject
{
private:
    // Track properties
    TLorentzVector TLV_part; // Lorentz vector properties of this particle

public:
    PythiaParticle() :
        TLV_part()
    {
    }
	~PythiaParticle()
	{
	}

	// setters
	void set_TLV_part(TLorentzVector tlv)     { TLV_part = tlv; }

	// getters
	TLorentzVector get_TLV_part() const       { return TLV_part;   }

	ClassDef(PythiaParticle,1);  // A simple track of a particle
};
//----------------------------------------------------------------------------------------





//----------------------------------------------------------------------------------------
class PythiaEvent : public TObject
{
private:
    Int_t eventNumber;
    Int_t   id; // id

    UShort_t      fNumParticles; // number of tracks in event

    TClonesArray* fParticles;      //->

public:
    PythiaEvent() :
        eventNumber(-1),id(-1),fNumParticles(0)
    {
        fParticles         = new TClonesArray( "PythiaParticle", 10 );
    }
	~PythiaEvent()
	{
	    delete fParticles;
            fParticles = NULL;
	}

        void       setEventNumber(Int_t e)          { eventNumber = e;                         }
	Int_t    getEventNumber() const             { return eventNumber;                      }

	void       setid(Int_t  r)                    { id = r;                        }
	Int_t      getid() const                      { return id;                     }

        //----------------------------
	PythiaParticle* createParticle()
	{
	    if (fNumParticles == fParticles->GetSize())
		fParticles->Expand( fNumParticles + 10 );
	    if (fNumParticles >= 50000)
	    {
		Fatal( "PythiaEvent::createParticle()", "ERROR: Too many tracks (>10000)!" );
		exit( 2 );
	    }

	    new((*fParticles)[fNumParticles++]) PythiaParticle;
	    return (PythiaParticle*)((*fParticles)[fNumParticles - 1]);
	}
	void clearParticleList()
	{
	    fNumParticles   = 0;
	    fParticles      ->Clear();
	}
	UShort_t getNumParticles() const
	{
	    return fNumParticles;
	}
	PythiaParticle* getParticle(UShort_t i) const
	{
	    return i < fNumParticles ? (PythiaParticle*)((*fParticles)[i]) : NULL;
        }
        //----------------------------



ClassDef(PythiaEvent,1);  // A simple event compiled of tracks
};
//----------------------------------------------------------------------------------------



#endif // __PYTHIAEVENT_H__