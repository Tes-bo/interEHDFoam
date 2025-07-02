if (cutted == false)
{
    if (equationType=="NP")
    {
        if(averageMethod=="linear")
        {
            if(cutEps)
            {
                if (limitFlux)
                {
                    epsSuf=(alphaSuf+fvc::interpolate(alpha1Sub)) * eps1 + (1 - (alphaSuf+fvc::interpolate(alpha1Sub))) * eps2;
                }
                else
                {
                    epsSuf = alphaSuf * eps1 + (1 - alphaSuf) * eps2;
                }           
            }
            else
            {
                epsSuf=fvc::interpolate(eps);
            }

            KSuf = alphaSuf * K1 + (1 - alphaSuf) * K2;
            DSuf = alphaSuf * D1 + (1 - alphaSuf) * D2;
        }
        else
        {
            if(cutEps)
            {            
                if (limitFlux)
                {
                    epsSuf = (eps1 * eps2) / (((alphaSuf+fvc::interpolate(alpha1Sub)) * eps2) + ((1 - (alphaSuf+fvc::interpolate(alpha1Sub))) * eps1));
                }
                else
                {
                    epsSuf = (eps1 * eps2) / ((alphaSuf * eps2) + ((1 - alphaSuf) * eps1));
                }
            }
            else
            {
                epsSuf=fvc::interpolate(eps);
            }
            KSuf = (K1 * K2) / ((alphaSuf * K2) + ((1 - alphaSuf) * K1));     
            DSuf = (D1 * D2) / ((alphaSuf * D2) + ((1 - alphaSuf) * D1)); 
        }
    }
    else
    {
        if(averageMethod=="linear")
        {
            if(cutEps)
            {
                if (limitFlux)
                {
                    epsSuf=(alphaSuf+fvc::interpolate(alpha1Sub)) * eps1 + (1 - (alphaSuf+fvc::interpolate(alpha1Sub))) * eps2;
                }
                else
                {
                    epsSuf = alphaSuf * eps1 + (1 - alphaSuf) * eps2;
                }           
            }
            else
            {
                epsSuf=fvc::interpolate(eps);
            }
            sigSuf = alphaSuf * sig1 + (1 - alphaSuf) * sig2;
        }
        else
        {
            if(cutEps)
            {            
                if (limitFlux)
                {
                    epsSuf = (eps1 * eps2) / (((alphaSuf+fvc::interpolate(alpha1Sub)) * eps2) + ((1 - (alphaSuf+fvc::interpolate(alpha1Sub))) * eps1));
                }
                else
                {
                    epsSuf = (eps1 * eps2) / ((alphaSuf * eps2) + ((1 - alphaSuf) * eps1));
                }
            }
            else
            {
                epsSuf=fvc::interpolate(eps);
            }
            sigSuf = (sig1 * sig2) / ((alphaSuf * sig2) + ((1 - alphaSuf) * sig1));
        }
    }

cutted =true;
}

