typedef float Real;

class parameterSet
{
public:
    parameterSet();

    Real getValue(String str)
    {
        return parameterMap[str];
    }

protected:
    virtual void instantiateValuesFrom(parameterSet*, parameterSet*);
    std::map<String, Real> parameterMap;
};


spatialParticleParameterSet()
{
public:
    enum KDExtrap {MIN, MAX, AV};

    spatialParticleParameterSet()
    {
        parameterMap["k_D"] = 0.0f;
        parameterMap["radius"] = 0.0f;
        parameterMap["mass"] = 0.0f;
    }

    void instantiateValuesFrom( parameterSet* ps1, parameterSet* ps2)
    {
        parameterMap["mass"] = (*massExtrap)(ps1->getValue("mass"),
                                             ps2->getValue("mass"));

        parameterMap["radius"] = (*radiusExtrap)(ps1->getValue("radius"),
                                                 ps2->getValue("radius"));

        parameterMap["k_D"] = (*k_DExtrap)(ps1->getValue("k_D"),
                                           ps2->getValue("k_D"));
    }

    static void setKDEExtrapolation(KDExtrap extrapCode)
    {
        if (extrapCode == MIN )
        {
            k_DExtrap = &spatialParticleParameterSet::extrapolateDiffusionCoeffMin;
        }
        else if (extrapCode == MAX )
        {
            k_DExtrap = &spatialParticleParameterSet::extrapolateDiffusionCoeffMax;
        }
        if (extrapCode == AV )
        {
            k_DExtrap = &spatialParticleParameterSet::extrapolateDiffusionCoeffAv;
        }

    }

    static void set(K)
    {
        
    }

    void setKDExtrap();
    void setMassExtrap();

protected:

    typedef Real (*FunctionPtr)(Real, Real) ExtrapolationFunctionPtr;

    static ExtrapolationFunctionPtr massExtrap;
    static ExtrapolationFunctionPtr radiusExtrap;
    static ExtrapolationFunctionPtr k_DExtrap;
    
    Real mass;
    Real radius;
    Real diffusion_coeff;

    static Real extrapolateDiffusionCoeffMin( Real diff1, Real diff2)
    {
        diff1 < diff2? diff1 : diff2;
    }

    static Real extrapolateDiffusionCoeffAv( Real diff1, Real diff2)
    {
        return (diff1 + diff2 ) / 2.0f;
    }

    static Real extrapolateDiffusionCoeffMax( Real diff1, Real diff2)
    {
        diff1 > diff2? diff1 : diff2;
    }
        
};

class spatialReactionParameterSet()
{
public:
};



class nonspatialParticleParameterSet()
{
public:
};

class nonspatialReactionParameterSet()
{
public:


};
