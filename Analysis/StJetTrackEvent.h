#ifndef __STJETTRACKEVENT_H__
#define __STJETTRACKEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector2.h"
//#include "StarClassLibrary/StThreeVectorF.hh"

// A. Schmah 05/31/2013
// A. Schmah 11/11/2013
// A. Schmah 02.03/2021 added ALICE EMCal container

//----------------------------------------------------------------------------------------
class StMCParticle : public TObject
{
private:
    TVector3 TV3_particle_vertex; //
    TLorentzVector TLV_particle;
    Short_t MC_charge;
    Int_t PDGcode;
    Int_t index_particle;
    Int_t index_mother;
    Int_t N_daughters;
    Int_t arr_index_daughters[5];
    Int_t MClabel;

public:
    StMCParticle() :
      TV3_particle_vertex(),TLV_particle(),MC_charge(0),PDGcode(0),index_particle(0),index_mother(0),N_daughters(0),arr_index_daughters(),MClabel(0)
    {

    }
        ~StMCParticle()
        {

        }

        void       set_TV3_particle_vertex(TVector3 TV3_particle_vertex_in)           { TV3_particle_vertex = TV3_particle_vertex_in;    }
        TVector3   get_TV3_particle_vertex() const                                    { return TV3_particle_vertex;                      }

        void set_TLV_particle(TLorentzVector tlv)                                     { TLV_particle = tlv;                              }
        TLorentzVector get_TLV_particle() const                                       { return TLV_particle;                             }

        void       set_MC_charge(Short_t ch)                                          { MC_charge = ch;                                  }
        Short_t    get_MC_charge() const                                              { return MC_charge;                                }

        void       set_PDGcode(Int_t PDGcode_in)                                      { PDGcode = PDGcode_in;                            }
        Int_t      get_PDGcode() const                                                { return PDGcode;                                  }

        void       set_MClabel(Int_t label_in)                                        { MClabel = label_in;                              }
        Int_t      get_MClabel() const                                                { return MClabel;                                  }

        void       set_index_mother(Int_t index_mother_in)                            { index_mother = index_mother_in;                  }
        Int_t      get_index_mother() const                                           { return index_mother;                             }

        void       set_index_particle(Int_t index_particle_in)                        { index_particle = index_particle_in;              }
        Int_t      get_index_particle() const                                         { return index_particle;                           }

        void       set_N_daughters(Int_t N_daughters_in)                              { N_daughters = N_daughters_in;                    }
        Int_t      get_N_daughters() const                                            { return N_daughters;                              }

        void       set_arr_index_daughters(Int_t index, Int_t index_daughter_in)      { arr_index_daughters[index] = index_daughter_in;  }
        Int_t      get_arr_index_daughters(Int_t index) const                         { return arr_index_daughters[index];               }


        ClassDef(StMCParticle,3);
};
//----------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------
class StEMCal : public TObject
{
    // EMCal data container
private:
    // Track properties
    TLorentzVector TLV_EMCal; // Lorentz vector properties of this EMCal cluster
    Float_t  cluster_pos[3];
    Float_t  Dx_to_track;
    Float_t  Dz_to_track;
    UShort_t N_cells_in_cluster;
    Bool_t   isExotic;
    Bool_t   isPHOS;
    Bool_t   isEMCal;

public:
    StEMCal() :
        TLV_EMCal(),cluster_pos(),Dx_to_track(0),Dz_to_track(0),N_cells_in_cluster(0),isExotic(0),isPHOS(0),isEMCal(0)
    {

    }
        ~StEMCal()
        {

        }

        // setters
        void set_TLV_EMCal(TLorentzVector tlv)                {TLV_EMCal  = tlv;    }
        void set_cluster_pos(Float_t x, Float_t y, Float_t z) {cluster_pos[0] = x; cluster_pos[1] = y; cluster_pos[2] = z;}
        void set_Dx_to_track(Float_t f)                       {Dx_to_track = f;}
        void set_Dz_to_track(Float_t f)                       {Dz_to_track = f;}
        void set_N_cells_in_cluster(UShort_t s)               {N_cells_in_cluster = s;}
        void set_isExotic(Bool_t b)                           {isExotic = b;}
        void set_isPHOS(Bool_t b)                             {isPHOS = b;}
        void set_isEMCal(Bool_t b)                            {isEMCal = b;}

        // getters
        TLorentzVector get_TLV_EMCal() const                  { return TLV_EMCal;   }
        Float_t        get_cluster_pos(Int_t index) const     { return cluster_pos[index];}
        Float_t        get_Dx_to_track() const                { return Dx_to_track;}
        Float_t        get_Dz_to_track() const                { return Dz_to_track;}
        UShort_t       get_N_cells_in_cluster() const         { return N_cells_in_cluster;}
        Bool_t         get_isExotic() const                   { return isExotic;}
        Bool_t         get_isPHOS() const                     { return isPHOS;}
        Bool_t         get_isEMCal() const                    { return isEMCal;}


        ClassDef(StEMCal,1)  //
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class StJetTrackParticle : public TObject
{
private:
    // Track properties
    //TLorentzVector TLV_Particle_prim; // Lorentz vector properties of this particle (primary track)

    Short_t        dca_to_prim_xy; // distance of closest approach to mother particle decay vertex (or primary vertex)
    Short_t        dca_to_prim_z; // distance of closest approach to mother particle decay vertex (or primary vertex)
    //Float_t        Particle_m2;
    //Float_t        Particle_nSigmaPi;
    //Float_t        Particle_nSigmaK;
    //Float_t        Particle_nSigmaP;
    //Float_t        Particle_qp;

    UShort_t       status; // status of track: bit 0: ITS refit, bit1: TPC refit
    Short_t        TPCchi2; // TPC chi2
    UShort_t       TPCdEdx; // Energy loss information of TPC
    UShort_t       TOFsignal; // Time-of-flight
    UShort_t       Track_length; // length of track
    //UShort_t     Length_in_active_zone; //the length of the track in cm in "active volume" of the TPC
    UShort_t       N_TPC_cls;  //number of used TPC clusters
    //Short_t        Particle_px;
    //Short_t        Particle_py;
    //Short_t        Particle_pz;
    //Int_t          ITS_layer;
    ULong64_t      Particle_p;
    Int_t          MClabel;

    Short_t particle_charge;
    

public:
    StJetTrackParticle() :
        //TLV_Particle_prim(),
        dca_to_prim_xy(-1),dca_to_prim_z(-1),
        //Particle_m2(-1),Particle_nSigmaPi(-1),Particle_nSigmaK(-1),Particle_nSigmaP(-1),Particle_qp(-1),
        status(0),TPCchi2(-3),
        TPCdEdx(-3),TOFsignal(-3),Track_length(-3),N_TPC_cls(0),Particle_p(0),particle_charge(0),MClabel(0)//,Particle_px(0),Particle_py(0),Particle_pz(0),//,ITS_layer(0)
    {

    }
        ~StJetTrackParticle()
        {

        }

        // setters
        void set_dca_to_prim(Float_t f1, Float_t f2)        { dca_to_prim_xy = (Short_t)(f1*100); dca_to_prim_z = (Short_t)(f2*100);           }
        //void set_Particle_m2(Float_t f)                     { Particle_m2 = f;            }
        //void set_Particle_nSigmaPi(Float_t f)               { Particle_nSigmaPi = f;      }
        //void set_Particle_nSigmaK(Float_t f)                { Particle_nSigmaK = f;       }
        //void set_Particle_nSigmaP(Float_t f)                { Particle_nSigmaP = f;       }
        //void set_Particle_qp(Float_t f)                     { Particle_qp = f;            }
        //void set_TLV_Particle_prim(TLorentzVector tlv)      { TLV_Particle_prim = tlv;    }

        //void set_px_Particle(Float_t f)                     { Particle_px = (Short_t)(f*200);          }
        //void set_py_Particle(Float_t f)                     { Particle_py = (Short_t)(f*200);          }
        //void set_pz_Particle(Float_t f)                     { Particle_pz = (Short_t)(f*200);          }
        void set_px_Particle(Float_t f)                     { Long64_t Particle_px = (f*6000);
                                                              if(Particle_px < 0)  Particle_p |=(ULong64_t)1 << 20;
                                                              if(Particle_px < 0)  Particle_px = (-1)*Particle_px;
                                                              for(Int_t i_bit = 0; i_bit < 20; i_bit++)
                                                              {
                                                                  if((Particle_px & ((Long64_t)1 <<  i_bit)))  Particle_p |= (ULong64_t)1 << i_bit;
                                                              };                                       }

        void set_py_Particle(Float_t f)                     { Long64_t Particle_py = (f*6000);
                                                              if(Particle_py < 0)  Particle_p |=(ULong64_t)1 << 41;
                                                              if(Particle_py < 0)  Particle_py = (-1)*Particle_py;
                                                              for(Int_t i_bit = 21; i_bit < 41 ;i_bit++)
                                                              {
                                                                  if((Particle_py & ((Long64_t)1 <<  (i_bit-21))))  Particle_p |= (ULong64_t)1 << i_bit;
                                                              };                                       }

        void set_pz_Particle(Float_t f)                     { Long64_t Particle_pz = (f*6000);
                                                              if(Particle_pz < 0)  Particle_p |=(ULong64_t)1 << 62;
                                                              if(Particle_pz < 0)  Particle_pz = (-1)*Particle_pz;
                                                              for(Int_t i_bit = 42; i_bit < 62 ;i_bit++)
                                                              {
                                                                  if((Particle_pz & ((Long64_t)1 <<  (i_bit-42))))  Particle_p |= (ULong64_t)1 << i_bit;
                                                              };                                       }

        void set_charge_Particle(Short_t ch)                { particle_charge = ch;                    }

        void  set_MClabel(Int_t label_in)                   { MClabel = label_in;                      }
        Int_t get_MClabel() const                           { return MClabel;                          }
       

        void setStatus(UShort_t s)                          { status = s;                              }
	void setTPCchi2(Float_t f)                          { TPCchi2 = (Short_t)(f*10);               }
	void setTPCdEdx(Float_t f)                          { TPCdEdx = (UShort_t)(f*10);              }
	void setTOFsignal(Float_t f)                        { TOFsignal = (UShort_t)(f/5);             }
        void setTrack_length(Float_t f)                     { Track_length = (UShort_t)(f*10);         }
        void setN_TPC_cls(UShort_t s)                       { N_TPC_cls = s;                           }

        // getters
        Float_t get_dca_to_prim_xy()             const        { return ((Float_t)dca_to_prim_xy)/100.0;        }
        Float_t get_dca_to_prim_z()              const        { return ((Float_t)dca_to_prim_z)/100.0;         }
        //Float_t get_Particle_m2 ()             const        { return Particle_m2;         }
        //Float_t get_Particle_nSigmaPi()        const        { return Particle_nSigmaPi;   }
        //Float_t get_Particle_nSigmaK()         const        { return Particle_nSigmaK;    }
        //Float_t get_Particle_nSigmaP()         const        { return Particle_nSigmaP;    }
        //Float_t get_Particle_qp()              const        { return Particle_qp;         }
        //TLorentzVector get_TLV_Particle_prim() const        { return TLV_Particle_prim.SetPxPyPzE(((Float_t)Particle_px)/200.0,((Float_t)Particle_py)/200.0,((Float_t)Particle_pz)/200.0,0.1349766);   }

        //Float_t get_px_Particle()              const        { return ((Float_t)Particle_px)/200.0;         }
        //Float_t get_py_Particle()              const        { return ((Float_t)Particle_py)/200.0;         }
        //Float_t get_pz_Particle()              const        { return ((Float_t)Particle_pz)/200.0;         }

        Float_t get_px_Particle()                           { Long64_t Particle_px = 0;
                                                              for(Int_t i_bit = 0; i_bit < 20 ;i_bit++)
                                                              {
                                                                  if((Particle_p & ((ULong64_t) 1 <<  i_bit))) Particle_px |= (Long64_t)1 << i_bit;
                                                              }
                                                              if((Particle_p & ((ULong64_t) 1 <<  20)))  Particle_px = (-1)*Particle_px;
                                                              return ((Float_t)Particle_px)/6000.0;                   }

        Float_t get_py_Particle()                           { Long64_t Particle_py = 0;
                                                              for(Int_t i_bit = 0; i_bit < 20 ;i_bit++)
                                                              {
                                                                  if((Particle_p & ((ULong64_t) 1 <<  (i_bit+21)))) Particle_py |= (Long64_t)1 << i_bit;
                                                              }
                                                              if((Particle_p & ((ULong64_t) 1 <<  41)))  Particle_py = (-1)*Particle_py;
                                                              return ((Float_t)Particle_py)/6000.0;                   }

        Float_t get_pz_Particle()                           { Long64_t Particle_pz = 0;
                                                              for(Int_t i_bit = 0; i_bit < 20 ;i_bit++)
                                                              {
                                                                  if((Particle_p & ((ULong64_t) 1 <<  (i_bit+42)))) Particle_pz |= (Long64_t)1 << i_bit;
                                                              }
                                                              if((Particle_p & ((ULong64_t) 1 <<  62)))  Particle_pz = (-1)*Particle_pz;
                                                              return ((Float_t)Particle_pz)/6000.0;                   }

        Short_t get_charge_Particle()          const        { return particle_charge;                      }
        


        UShort_t getStatus()                   const        { return status;                               }
        Float_t  getTPCchi2()                  const        { return ((Float_t)TPCchi2)/10.0;              }
	Float_t  getTPCdEdx()                  const        { return ((Float_t)TPCdEdx)/10.0;              }
	Float_t  getTOFsignal()                const        { return ((Float_t)TOFsignal)*5.0;             }
        Float_t  getTrack_length()             const        { return ((Float_t)Track_length)/10.0;         }
        UShort_t getN_TPC_cls()                const        { return N_TPC_cls;                            }

        Int_t    ITSLayer()                                 { Int_t ITS_layer = 0;
                                                              for(Int_t i_ITS_layer = 0; i_ITS_layer < 6; ++i_ITS_layer)
                                                              {
                                                                  if((status & ((UShort_t)1 << i_ITS_layer)) )  ITS_layer+=1;
                                                              }
                                                              return ITS_layer;                                       }

        Bool_t   IsITSrefit()                               { if((status & ((UShort_t)1 <<  6))) return kTRUE; else return kFALSE; }
        Bool_t   IsBit32()                                  { if((status & ((UShort_t)1 <<  7))) return kTRUE; else return kFALSE; }
        Bool_t   IsBit64()                                  { if((status & ((UShort_t)1 <<  8))) return kTRUE; else return kFALSE; }
        Bool_t   IsBit128()                                 { if((status & ((UShort_t)1 <<  9))) return kTRUE; else return kFALSE; }
        Bool_t   IsBit256()                                 { if((status & ((UShort_t)1 << 10))) return kTRUE; else return kFALSE; }
        Bool_t   IsBit512()                                 { if((status & ((UShort_t)1 << 11))) return kTRUE; else return kFALSE; }

        ClassDef(StJetTrackParticle,3)  // A simple track of a particle
};
//------------------------------------------------------------------------------------



//------------------------------------------------------------------------------------
class StJetTrackEvent : public TObject
{
private:
    Float_t x;
    Float_t y;
    Float_t z;
    Int_t   id;
    Float_t mult;
    Int_t   n_prim;
    Int_t   n_TRD_tracklets;
    Int_t   n_tof_hits;
    Int_t   SE_ME_flag;
    Int_t   n_pileup_vertices_tracks;
    Float_t multiplicity_V0A;
    Float_t multiplicity_V0C;
    Float_t multiplicity_V0A_cell[32];
    Float_t multiplicity_V0C_cell[32];
    Short_t flag_Data_MC; // 0 = data, 1 = MC
    Int_t   mult_used_for_ME;
    Float_t Psi_EP;
  //ULong64_t trigger_mask;
  Bool_t is_kINT7;
  Bool_t is_kCentral;
  Bool_t is_kSemiCentral;
  Float_t Psi2_EP_V0;
  Float_t Psi3_EP_V0;

    TVector2 QvecEtaPos;
    TVector2 QvecEtaNeg;

    Float_t   cent_class_ZNA; // ZDC neutral A
    Float_t   cent_class_ZNC; // ZDC neutral C
    Float_t   cent_class_V0A; // V0 A
    Float_t   cent_class_V0C; // V0 C
    Float_t   cent_class_V0M; // V0 average
    Float_t   cent_class_CL0; // clusters in layer 0
    Float_t   cent_class_CL1; // clusters in layer 1
    Float_t   cent_class_SPD; // SPD
    Float_t   cent_class_V0MEq; //
    Float_t   cent_class_V0AEq; //
    Float_t   cent_class_V0CEq; //

    Float_t BeamIntAA; // ZDC coincidence rate
    Float_t T0zVertex; // z-vertex position from VPD

    TString TriggerWord; // Trigger word

    UShort_t      fNumParticle;
    UShort_t      fNumEMCal;
    Int_t         fNumMCparticles; // number of MC particles

    TClonesArray* fParticle; //->
    TClonesArray* fMCparticles;  //->
    TClonesArray* fEMCal; //->

public:
    StJetTrackEvent() :
        x(-1),y(-1),z(-1),id(-1),mult(0),n_prim(0),n_TRD_tracklets(0),
        n_tof_hits(0),n_pileup_vertices_tracks(0),multiplicity_V0A(0),multiplicity_V0C(0),multiplicity_V0A_cell(),multiplicity_V0C_cell(),SE_ME_flag(-1),QvecEtaPos(),QvecEtaNeg(),
        cent_class_ZNA(0),cent_class_ZNC(0),cent_class_V0A(0),cent_class_V0C(0),cent_class_V0M(0),cent_class_CL0(0),cent_class_CL1(0),
        cent_class_SPD(0),cent_class_V0MEq(0),cent_class_V0AEq(0),cent_class_V0CEq(0),BeamIntAA(-1),T0zVertex(-1),flag_Data_MC(0),
        TriggerWord(),fNumParticle(0),fNumEMCal(0),fNumMCparticles(0),fMCparticles(),mult_used_for_ME(0),Psi_EP(0),is_kINT7(0),is_kCentral(0),is_kSemiCentral(0),Psi2_EP_V0(9999),Psi3_EP_V0(9999)
    {
        fParticle      = new TClonesArray( "StJetTrackParticle", 10 );
        fEMCal         = new TClonesArray( "StEMCal", 10 );
        fMCparticles    = new TClonesArray( "StMCParticle", 10 );
    }
        ~StJetTrackEvent()
        {
            delete fParticle;
            fParticle = NULL;
            delete fMCparticles;
            fMCparticles = NULL;
        }

  //void      set_trigger_mask(ULong64_t r)       { trigger_mask = r;              }
  //ULong64_t get_trigger_mask() const            { return trigger_mask;           }
        void       set_is_kINT7(Bool_t r)             { is_kINT7 = r;                  }
        Bool_t     get_is_kINT7() const               { return is_kINT7;               }

        void       set_is_kCentral(Bool_t r)          { is_kCentral = r;               }
        Bool_t     get_is_kCentral() const            { return is_kCentral;            }

        void       set_is_kSemiCentral(Bool_t r)      { is_kCentral = r;               }
        Bool_t     get_is_kSemiCentral() const        { return is_kCentral;            }

        void       setx(Float_t r)                    { x = r;                         }
        Float_t    getx() const                       { return x;                      }

        void       sety(Float_t r)                    { y = r;                         }
        Float_t    gety() const                       { return y;                      }

        void       setz(Float_t r)                    { z = r;                         }
        Float_t    getz() const                       { return z;                      }

        void       setid(Int_t  r)                    { id = r;                        }
        Int_t      getid() const                      { return id;                     }

        void       setmult(Float_t r)                 { mult = r;                      }
        Float_t    getmult() const                    { return mult;                   }

        void       setn_prim(Int_t r)                 { n_prim = r;                    }
        Int_t      getn_prim() const                  { return n_prim;                 }

        void       setn_TRD_tracklets(Int_t r)        { n_TRD_tracklets = r;           }
        Int_t      getn_TRD_tracklets() const         { return n_TRD_tracklets;        }

        void       setn_tof_hits(Int_t r)             { n_tof_hits = r;                }
        Int_t      getn_tof_hits() const              { return n_tof_hits;             }

        void       setn_pileup_vertices_tracks(Int_t r) { n_pileup_vertices_tracks = r;                }
        Int_t      getn_pileup_vertices_tracks() const  { return n_pileup_vertices_tracks;             }

        void       setmult_V0A(Float_t r)             { multiplicity_V0A = r;          }
        Float_t    getmult_V0A() const                { return multiplicity_V0A;       }

        void       setmult_V0C(Float_t r)             { multiplicity_V0C = r;          }
        Float_t    getmult_V0C() const                { return multiplicity_V0C;       }

        void       setmult_V0A_cell(Int_t ich, Float_t r)   { multiplicity_V0A_cell[ich] = r;    }
        Float_t    getmult_V0A_cell(Int_t ich) const        { return multiplicity_V0A_cell[ich]; }

        void       setmult_V0C_cell(Int_t ich, Float_t r)   { multiplicity_V0C_cell[ich] = r;    }
        Float_t    getmult_V0C_cell(Int_t ich) const        { return multiplicity_V0C_cell[ich]; }

        void       setmult_used_for_ME(Int_t  r)      { mult_used_for_ME = r;          }
        Int_t      getmult_used_for_ME() const        { return mult_used_for_ME;       }

        void       setPsi_EP(Float_t r)               { Psi_EP = r;                    }
        Float_t    getPsi_EP() const                  { return Psi_EP;                 }

        void       setPsi2_EP_V0(Float_t r)           { Psi2_EP_V0 = r;                    }
        Float_t    getPsi2_EP_V0() const              { return Psi2_EP_V0;                 }

        void       setPsi3_EP_V0(Float_t r)           { Psi3_EP_V0 = r;                    }
        Float_t    getPsi3_EP_V0() const              { return Psi3_EP_V0;                 }

        void       setSE_ME_flag(Int_t r)             { SE_ME_flag = r;                }
        Int_t      getSE_ME_flag() const              { return SE_ME_flag;             }

        void       setQvecEtaPos(TVector2 vec)        { QvecEtaPos = vec;              }
        TVector2   getQvecEtaPos() const              { return QvecEtaPos;             }

        void       setQvecEtaNeg(TVector2 vec)        { QvecEtaNeg = vec;              }
        TVector2   getQvecEtaNeg() const              { return QvecEtaNeg;             }

        void       setcent_class_ZNA(Float_t r)             { cent_class_ZNA = r;                }
	Float_t      getcent_class_ZNA() const              { return cent_class_ZNA;             }

	void       setcent_class_ZNC(Float_t r)             { cent_class_ZNC = r;                }
	Float_t      getcent_class_ZNC() const              { return cent_class_ZNC;             }

	void       setcent_class_V0A(Float_t r)             { cent_class_V0A = r;                }
	Float_t      getcent_class_V0A() const              { return cent_class_V0A;             }

	void       setcent_class_V0C(Float_t r)             { cent_class_V0C = r;                }
	Float_t      getcent_class_V0C() const              { return cent_class_V0C;             }

	void       setcent_class_V0M(Float_t r)             { cent_class_V0M = r;                }
	Float_t      getcent_class_V0M() const              { return cent_class_V0M;             }

	void       setcent_class_CL0(Float_t r)             { cent_class_CL0 = r;                }
	Float_t      getcent_class_CL0() const              { return cent_class_CL0;             }

	void       setcent_class_CL1(Float_t r)             { cent_class_CL1 = r;                }
	Float_t      getcent_class_CL1() const              { return cent_class_CL1;             }

	void       setcent_class_SPD(Float_t r)             { cent_class_SPD = r;                }
	Float_t      getcent_class_SPD() const              { return cent_class_SPD;             }

	void       setcent_class_V0MEq(Float_t r)             { cent_class_V0MEq = r;                }
	Float_t      getcent_class_V0MEq() const              { return cent_class_V0MEq;             }

	void       setcent_class_V0AEq(Float_t r)             { cent_class_V0AEq = r;                }
	Float_t      getcent_class_V0AEq() const              { return cent_class_V0AEq;             }

	void       setcent_class_V0CEq(Float_t r)             { cent_class_V0CEq = r;                }
	Float_t      getcent_class_V0CEq() const              { return cent_class_V0CEq;             }

	void       setBeamIntAA(Float_t r)                 { BeamIntAA = r;                      }
	Float_t    getBeamIntAA() const                    { return BeamIntAA;                   }

	void       setT0zVertex(Float_t r)            { T0zVertex = r;                     }
	Float_t    getT0zVertex() const               { return T0zVertex;                  }

        void       setTriggerWord(TString s)          { TriggerWord = s;}
        TString    getTriggerWord() const             { return TriggerWord; }

        void       set_flag_Data_MC(UShort_t s)             { flag_Data_MC = s;                }
	UShort_t   get_flag_Data_MC() const                 { return flag_Data_MC;             }

        //--------------------------------------
        // Particle
        StJetTrackParticle* createParticle()
        {
            if (fNumParticle == fParticle->GetSize())
                fParticle->Expand( fNumParticle + 5 );
            if (fNumParticle >= 50000)
            {
                Fatal( "StJetTrackParticle::createParticle()", "ERROR: Too many Particles (>5000)!" );
                exit( 2 );
            }

            new((*fParticle)[fNumParticle++]) StJetTrackParticle;
            return (StJetTrackParticle*)((*fParticle)[fNumParticle - 1]);
        }
        void clearParticleList()
        {
            fNumParticle   = 0;
            fParticle      ->Clear();
        }
        UShort_t getNumParticle() const
        {
            return fNumParticle;
        }
        StJetTrackParticle* getParticle(UShort_t i) const
        {
            return i < fNumParticle ? (StJetTrackParticle*)((*fParticle)[i]) : NULL;
        }
        //--------------------------------------



        //--------------------------------------
        // EMCal
        StEMCal* createEMCal()
        {
            if (fNumEMCal == fEMCal->GetSize())
                fEMCal->Expand( fNumEMCal + 5 );
            if (fNumEMCal >= 50000)
            {
                Fatal( "StEMCal::createEMCal()", "ERROR: Too many EMCals (>5000)!" );
                exit( 2 );
            }

            new((*fEMCal)[fNumEMCal++]) StEMCal;
            return (StEMCal*)((*fEMCal)[fNumEMCal - 1]);
        }
        void clearEMCalList()
        {
            fNumEMCal   = 0;
            fEMCal      ->Clear();
        }
        UShort_t getNumEMCal() const
        {
            return fNumEMCal;
        }
        StEMCal* getEMCal(UShort_t i) const
        {
            return i < fNumEMCal ? (StEMCal*)((*fEMCal)[i]) : NULL;
        }
        //----------------------------



        //----------------------------
        // Monte Carlo particle
        StMCParticle* createMCparticle() 
	{
            if (fNumMCparticles == fMCparticles->GetSize())
		fMCparticles->Expand( fNumMCparticles + 10 );
	    if (fNumMCparticles >= 650000)
	    {
		Fatal( "StMCParticle::createMCparticle()", "ERROR: Too many MC particles (>650000)!" );
		exit( 2 );
	    }

	    new((*fMCparticles)[fNumMCparticles++]) StMCParticle;
	    return (StMCParticle*)((*fMCparticles)[fNumMCparticles - 1]);
	}
	void clearMCparticleList()
	{
	    fNumMCparticles   = 0;
	    fMCparticles      ->Clear();
	}
	Int_t getMCparticles() const
	{
	    return fNumMCparticles;
	}
	StMCParticle* getMCparticle(Int_t i) const
	{
	    return i < fNumMCparticles ? (StMCParticle*)((*fMCparticles)[i]) : NULL;
        }
        //--------------------------------------


        ClassDef(StJetTrackEvent,2)  // A simple event compiled of tracks
};
//------------------------------------------------------------------------------------



#endif // __STJETTRACKEVENT_H__
