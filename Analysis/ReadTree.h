
Float_t get_px_Particle(Long64_t Particle_p)
{
    Long64_t Particle_px = 0;
    for(Int_t i_bit = 0; i_bit < 20 ;i_bit++)
    {
        if((Particle_p & ((ULong64_t) 1 <<  i_bit))) Particle_px |= (Long64_t)1 << i_bit;
    }
    if((Particle_p & ((ULong64_t) 1 <<  20)))  Particle_px = (-1)*Particle_px;
    return ((Float_t)Particle_px)/6000.0;
}

Float_t get_py_Particle(Long64_t Particle_p)
{
    Long64_t Particle_py = 0;
    for(Int_t i_bit = 0; i_bit < 20 ;i_bit++)
    {
        if((Particle_p & ((ULong64_t) 1 <<  (i_bit+21)))) Particle_py |= (Long64_t)1 << i_bit;
    }
    if((Particle_p & ((ULong64_t) 1 <<  41)))  Particle_py = (-1)*Particle_py;
    return ((Float_t)Particle_py)/6000.0;
}

Float_t get_pz_Particle(Long64_t Particle_p)
{
    Long64_t Particle_pz = 0;
    for(Int_t i_bit = 0; i_bit < 20 ;i_bit++)
    {
        if((Particle_p & ((ULong64_t) 1 <<  (i_bit+42)))) Particle_pz |= (Long64_t)1 << i_bit;
    }
    if((Particle_p & ((ULong64_t) 1 <<  62)))  Particle_pz = (-1)*Particle_pz;
    return ((Float_t)Particle_pz)/6000.0;
}

