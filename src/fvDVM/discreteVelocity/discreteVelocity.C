/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvDVM.H"
#include "discreteVelocity.H"
#include "fixedGradientFvPatchField.H"

#include "farFieldFvPatchField.H"
#include "farFieldFvsPatchField.H"

#include "bounceBackFvPatchField.H"
#include "bounceBackFvsPatchField.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::discreteVelocity::discreteVelocity
(
    fvDVM& dvm,
    const fvMesh& mesh,
    const Time& time,
    const scalar weight,
    const dimensionedVector xi,
    const label DVid,
    const label symXtargetDVid,
    const label symYtargetDVid,
    const label symZtargetDVid
)
:
    dvm_(dvm),
    mesh_(mesh),
    time_(time),
    weight_(weight),
    xi_(xi),
    myDVid_(DVid),
    bDVid_(dvm_.nXi() - DVid - 1),
    symXtargetDVid_(symXtargetDVid),
    symYtargetDVid_(symYtargetDVid),
    symZtargetDVid_(symZtargetDVid),
    gVol_
    (
        IOobject
        (
            "gVol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass/pow3(dimLength), 0.0
        )
    ),
    //BKG to DUGKS : add a field of supporting distribution function
    gBarVol_ 
    (
    	IOobject
    	(
    		"gBarVol" + name(DVid),
    		mesh_.time().timeName(),
    		mesh_,
    		IOobject::NO_READ,
    		IOobject::NO_WRITE
    	),
    	mesh_,
    	dimensionedScalar
        (
            "0", dimMass/pow3(dimLength), 0.0
        ),
        "fixedGradient"
    ),    
    gSurf_
    (
        IOobject
        (
            "gSurf" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass/pow3(dimLength), 0.0
        )
    ),
    gBarGrad_
    (
        IOobject
        (
            "gBarGrad" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "0", gVol_.dimensions()/dimLength, vector(0,0,0)
        ),
        "zeroGradient"
    )
{
    initDFtoEq();
    setBCtype();
    initBoundaryField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::discreteVelocity::~discreteVelocity()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::discreteVelocity::initDFtoEq()
{
    equilibrium
    (
        gVol_,
        dvm_.rhoVol(), 
        dvm_.Uvol() 
    );
}

void Foam::discreteVelocity::setBCtype()
{
    // Only from rho's B.C can we determine surfaceScalarField f's B.C.
    // for all patchi of gbarP, set to zeroGradient
    // May improve for inlet ingoing DF

     const GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = dvm_.rhoVol().boundaryField();

    forAll(rhoBCs, patchi)
    {
        if (rhoBCs[patchi].type() == "fixedValue") 
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "mixed", mesh_.boundary()[patchi], gSurf_
                )
            );
        }
        else if (rhoBCs[patchi].type() == "farField") 
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "farField", mesh_.boundary()[patchi], gSurf_
                )
            );
        }
        else if (rhoBCs[patchi].type() == "zeroGradient") 
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "zeroGradient", mesh_.boundary()[patchi], gSurf_
                )
            );
        }
        else if (rhoBCs[patchi].type() == "bounceBack") 
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "bounceBack", mesh_.boundary()[patchi], gSurf_
                )
            );
        }
        else if (rhoBCs[patchi].type() == "pressureIn") 
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "farField", mesh_.boundary()[patchi], gSurf_
                )
            );
        }
        else if (rhoBCs[patchi].type() == "pressureOut") 
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "farField", mesh_.boundary()[patchi], gSurf_
                )
            );
        }
    }
}

void Foam::discreteVelocity::initBoundaryField()
{
    GeometricField<scalar, fvsPatchField, surfaceMesh>::GeometricBoundaryField& 
        gBCs = gSurf_.boundaryField();

    forAll(gBCs, patchi)
    {
        if(gBCs[patchi].type() == "mixed" )
        {
            equilibrium
            (
                gBCs[patchi],
                dvm_.rhoVol().boundaryField()[patchi],
                dvm_.Uvol().boundaryField()[patchi]
            );
        }
    }
}

void Foam::discreteVelocity::updateGbarvol()
{
    //- get delta t
    static volScalarField gEq
    (
        IOobject
        (
            "gEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", gVol_.dimensions(), 0)
    );

    //- get gEq
    equilibrium
    (
        gEq,
        dvm_.rhoVol(), 
        dvm_.Uvol() 
    );
    gBarVol_ = (1.0 - dvm_.omega1())*gVol_ + dvm_.omega1()*gEq;
    gBarVol_.correctBoundaryConditions(); 
}

void Foam::discreteVelocity::updateGbarGrad()
{
      // Pre-setup
    // 1. correct the boundary value of gBarPvol_
    //    the gradient at boundary is known
    // 2. get the gradient
    //
    gBarGrad_ = fvc::grad(gBarVol_); 
    // The DVMsymmetry is  rocessed automatically in fvc::grad operator

    //
    // 3. correct the boundary value of the grad field
    //    to be used at next time
    gBarGrad_.correctBoundaryConditions();

    // 4. patch the normal component of boundary value of the grad 
    //    to the gradient field of the fixed 
    //    gradient feild of the gBarPvol_ 
    //    NOTE: we need the surfaceNormal Gradient
}

void Foam::discreteVelocity::updateGbarsurf()
{
    forAll(gBarGrad_.boundaryField(), patchi)
    {
        const vectorField n
        (
            mesh_.Sf().boundaryField()[patchi]
           /mesh_.magSf().boundaryField()[patchi]
        );

        if ( 
               gBarVol_.boundaryField()[patchi].type() != "empty" 
            && gBarVol_.boundaryField()[patchi].type() != "processor"
            && gBarVol_.boundaryField()[patchi].type() != "processorCyclic"
            && gBarVol_.boundaryField()[patchi].type() != "cyclic"
            && gBarVol_.boundaryField()[patchi].type() != "symmetryPlane"
           ) // only for fixed gradient gBarPvol
        {
            // normal component of the grad field
            fixedGradientFvPatchField<scalar>& gBarVolPatch = 
                refCast<fixedGradientFvPatchField<scalar> >
                (gBarVol_.boundaryField()[patchi]);

            forAll(gBarVolPatch, pFacei)
            {
                gBarVolPatch.gradient()[pFacei] =
                    gBarGrad_.boundaryField()[patchi][pFacei]&n[pFacei];
            }
        }
    }

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField& Cf = mesh_.Cf();
    const vectorField& Sf = mesh_.Sf();
    const vectorField& C = mesh_.C();

    const Field<scalar>& iGbarVol = gBarVol_;
    const Field<vector>& iGbarGrad = gBarGrad_;
    const Field<scalar>& ibGbarVol = dvm_.DVi(bDVid_).gBarVol();
    const Field<vector>& ibGbarGrad = dvm_.DVi(bDVid_).gBarGrad();

    // This is what we want to update in this function
    Field<scalar>& iGsurf = gSurf_;

    scalar dt = time_.deltaTValue(); 
    vector xii = xi_.value(); 
    label D = mesh_.nSolutionD();

    // internal faces first
    forAll(owner, facei)                    
    {
        label own = owner[facei];
        label nei = neighbour[facei];
        if ((xii&Sf[facei]) >= VSMALL) // comming from own
        {

            iGsurf[facei] = iGbarVol[own] 
              + (iGbarGrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));
        }
        else if ((xii&Sf[facei]) < -VSMALL) // comming form nei
        {
            iGsurf[facei] = iGbarVol[nei]
              + (iGbarGrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
        }
        else 
        {
            iGsurf[facei] = 0.5*
            (
                iGbarVol[nei] + ((iGbarGrad[nei])
               &(Cf[facei] - C[nei] - 0.5*xii*dt))
              + iGbarVol[own] + ((iGbarGrad[own])
               &(Cf[facei] - C[own] - 0.5*xii*dt))
            );
        }
    }

    // boundary faces
    forAll(gSurf_.boundaryField(), patchi)
    {
        word type = gSurf_.boundaryField()[patchi].type();
        fvsPatchField<scalar>& gSurfPatch = gSurf_.boundaryField()[patchi];
        const fvPatchField<vector>& Upatch = dvm_.Uvol().boundaryField()[patchi];
        const fvPatchField<scalar>& rhoPatch = dvm_.rhoVol().boundaryField()[patchi];
        const fvsPatchField<vector>& SfPatch =
            mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch =
            mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();
        
        //- NOTE: outging DF can be treate unifiedly for all BCs, including processor BC
        if (type == "zeroGradient")
        {
            gSurfPatch == gBarVol_.boundaryField()[patchi].patchInternalField();
        }
        else if (type == "mixed")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                //out or in ?
                if ((xii&SfPatch[facei]) > 0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarVol[faceCells[facei]]  
                      + ((iGbarGrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                }
                //incoming and parallel to face, not changed.    
            }
        }
        else if (type == "farField")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                //out or in ?
                if ((xii&SfPatch[facei]) > 0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarVol[faceCells[facei]] 
                      + ((iGbarGrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                }
                //due to update
                else // incomming, set to be equlibrium, give rho , U
                {
                    scalar CsSqr = dvm_.CsSqr().value();
                    scalar rho   = rhoPatch[facei];
                    vector U     = dvm.Uvol()[pOwner[facei]];

                    gSurfPatch[facei] =  weight_*rho*(1.0 + (xii&U)/CsSqr 
                                      + 0.5*((xii&U)*(xii&U))/CsSqr/CsSqr 
                                      - 0.5*magSqr(U)/CsSqr);
                }
            }
        }
        
        else if (type == "bounceBack")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                //out or in ?
                if((xii&SfPatch[facei]) < -VSMALL ) // incomming
                {
                    gSurfPatch[facei] = ibGbarVol[faceCells[facei]] 
                      + (ibGbarGrad[faceCells[facei]]
                       &(CfPatch[facei] - C[faceCells[facei]] + 0.5*xii*dt))
                     + 2.0*weight_*rhoPatch[facei]*(xii&Upatch[facei])/dvm_.CsSqr().value();//no slip wall;
                }
                else //if ((xii&SfPatch[facei]) > VSMALL ) // outgoing
                {
                    gSurfPatch[facei] = iGbarVol[faceCells[facei]] 
                      + ((iGbarGrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));

                //it's after collision,when outgoing, not changed.
                }
            }
        }

        else if ( type == "processor"
               || type == "cyclic"
               || type == "processorCyclic"
                 ) // parallel
        {
            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) > VSMALL) // outgoing
                {
                    gSurfPatch[facei] = iGbarVol[faceCells[facei]] 
                      + ((iGbarGrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt));
                } 
                else if ((xii&Sf[facei]) < -VSMALL) //incomming from processor boundaryField
                {
                    gSurfPatch[facei] = gBarVol_.boundaryField()[patchi][facei]
                      + ((gBarGrad_.boundaryField()[patchi][facei])
                       &(CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt));
                }
                else 
                {
                    gSurfPatch[facei] = 0.5*(
                            iGbarVol[faceCells[facei]] + 
                            (   (iGbarGrad[faceCells[facei]]) & (CfPatch[facei] - C[faceCells[facei]] - 0.5*xii*dt) )
                            + gBarVol_.boundaryField()[patchi][facei] + 
                            (   (gBarGrad_.boundaryField()[patchi][facei]) & (CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - 0.5*xii*dt) )
                    );
                }
            }
        }
    }
}
    
void Foam::discreteVelocity::updateGsurf()
{
    	
    static surfaceScalarField gSurfEq_
    (
        IOobject
        (
            "gSurfEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", gVol_.dimensions(), 0)
    );
    equilibrium
    (
        gSurfEq_,
        dvm_.rhoSurf(),
        dvm_.Usurf()      
    );
    //get the distribution function at internalface
    //Note that the boundary faces have not been calculated!
    gSurf_  = (1.0 - dvm_.omega2())* gSurf_ + dvm_.omega2()* gSurfEq_;

    //boundary field
    vector xii = xi_.value();
    // We do it mannuly for the outgoing DF!
    forAll(gSurf_.boundaryField(), patchi)
    {
        fvsPatchScalarField& gSurfPatch   = gSurf_.boundaryField()[patchi];
        fvsPatchScalarField& gSurfEqPatch = gSurfEq_.boundaryField()[patchi];

        const fvsPatchVectorField& SfPatch = mesh_.Sf().boundaryField()[patchi];

        //NOTE:  If keep the below block, for the Sod problem, -parallel running would be different from serial run.
        //       I also doubt that for cyclic boundary, the same problem will happen.
        //       But if I totally deleted it, my experience shows for hypersonic cylinder flow, the free stream BC will have problem.
        //       So, I keep it, but exclude the cases of processor/cyclic/processorCyclic.
        //       The reason I guess is that, for those 'coupled' type boundary condition, the boundary field will be calculated when 
        //       doing GeometricField calculation, such as the relaxiation computation. While for those boundary condition I've defined, 
        //       such as mixed or maxwell, I have to explicitly calculate the unchanged boundary values.
        forAll(gSurfPatch, facei)
        {
            if ( (xii&(SfPatch[facei])) > 0    // Here, if delted , the free stream BC may bo s problem
            && gBarVol_.boundaryField()[patchi].type() != "processor"
            && gBarVol_.boundaryField()[patchi].type() != "processorCyclic"
            && gBarVol_.boundaryField()[patchi].type() != "cyclic"
                  )
            //&& gSurf_.boundaryField()[patchi].type() != "DVMsymmetry")
            {
                gSurfPatch  = (1.0 - dvm_.omega2())* gSurfPatch + dvm_.omega2()* gSurfEqPatch;
            }
        }
    }
}

void Foam::discreteVelocity::updateGnewVol()
{
    // store the gTildePlus in g
    gVol_ = -1.0/3*gVol_ + 4.0/3*gBarVol_;

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField Sf = mesh_.Sf();
    const scalarField V = mesh_.V();
    const scalar dt = time_.deltaTValue();
    const vector xii = xi_.value();

    // internal faces
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        gVol_[own] -= ((xii&Sf[facei])*gSurf_[facei]*dt/V[own]);
        gVol_[nei] += ((xii&Sf[facei])*gSurf_[facei]*dt/V[nei]);
    }
    // boundary faces
    forAll(gSurf_.boundaryField(), patchi)
    {
        const fvsPatchField<scalar>& gSurfPatch =
            gSurf_.boundaryField()[patchi];
        const fvsPatchField<vector>& SfPatch =
            mesh_.Sf().boundaryField()[patchi];
        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
        forAll(gSurfPatch, pFacei)
        {
            const label own = pOwner[pFacei];
            gVol_[own] -= ((xii&SfPatch[pFacei]) *gSurfPatch[pFacei]*dt/V[own]);
        }
    }

}

void Foam::discreteVelocity::equilibrium
(
    volScalarField& geq,
    const volScalarField&  rho,
    const volVectorField&  U
)
{
    dimensionedVector xii = xi_;
    dimensionedScalar CsSqr = dvm_.CsSqr();
    geq =  weight_*rho*(1.0 + (xii&U)/CsSqr + 0.5*((xii&U)*(xii&U))/CsSqr/CsSqr - 0.50*magSqr(U)/CsSqr);
}

void Foam::discreteVelocity::equilibrium
(
    surfaceScalarField& geq,
    const surfaceScalarField&  rho,
    const surfaceVectorField&  U
)
{
    dimensionedVector xii = xi_;
    dimensionedScalar CsSqr = dvm_.CsSqr();
    geq =  weight_*rho*(1.0 + (xii&U)/CsSqr + 0.5*((xii&U)*(xii&U))/CsSqr/CsSqr - 0.50*magSqr(U)/CsSqr);
}

void Foam::discreteVelocity::equilibrium
(
    fvsPatchField<scalar>& geq,
    const fvPatchField<scalar>&  rho,
    const fvPatchField<vector>&  U
)
{
    vector xii = xi_.value();
    scalar CsSqr = dvm_.CsSqr().value();
    geq =  weight_*rho*(1.0 + (xii&U)/CsSqr + 0.5*((xii&U)*(xii&U))/CsSqr/CsSqr - 0.50*magSqr(U)/CsSqr);
}
// ************************************************************************* //
