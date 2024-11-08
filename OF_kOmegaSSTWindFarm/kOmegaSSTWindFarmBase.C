/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

    This is the openFoam kOmegaSST model corrected with Explicit Algebraic Reynolds Stress Correction Model (SecondaryFlow 2D)
    which is optimized for prediction of secondary flow in Square-Duct Reb = 3500.
    2023, A. Amarloo
    The default model coefficients  are
    \verbatim
        kOmegaSSTWindFarmBaseCoeffs
        {
            //original coefficients of KOSST
            alphaK1         0.85;
            alphaK2         1.0;
            alphaOmega1     0.5;
            alphaOmega2     0.856;
            beta1           0.075;
            beta2           0.0828;
            betaStar        0.09;
            gamma1          5/9;
            gamma2          0.44;
            a1              0.31;
            b1              1.0;
            c1              10.0;
            F3              no;


            //SecondaryFlow coefficients
            SecondaryFlow          ModelII; \\ off | ModelI | ModelII
            C0I             0.1;
            C1I             0.1;
            C2I             0.1;
            C0II            0.1;
            C1II            0.1;
            C2II            0.1;

            // Optional decay control
            decayControl    yes;
            kInf            \<far-field k value\>;
            omegaInf        \<far-field omega value\>;
        }
    \endverbatim



\*---------------------------------------------------------------------------*/

#include "kOmegaSSTWindFarmBase.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::F3() const
{
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


template<class BasicEddyViscosityModel>
void kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2
)
{
    // Correct the turbulence viscosity
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


template<class BasicEddyViscosityModel>
void kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicEddyViscosityModel>
tmp<fvVectorMatrix> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::divDevReff
(
    volVectorField& U
) const
{
    Info << "In: nonlinearRST::divDevReff()" << endl;
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)    // linear part
      + fvc::div(dev(2.*this->k_*this->bijDelta_))  // non-linear part
    );
}

template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    return betaStar_*omega_();
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omega_()
       *max(a1_*omega_(), b1_*F2*sqrt(S2))
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::kOmegaSSTWindFarmBase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0
        )
    ),
    SecondaryFlow_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "SecondaryFlow",
            this->coeffDict_,
            0
        )
    ),
    C0I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C0I",
            this->coeffDict_,
            0
        )
    ),
    C1I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C1I",
            this->coeffDict_,
            0
        )
    ),
    C2I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2I",
            this->coeffDict_,
            0
        )
    ),
    Separation_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Separation",
            this->coeffDict_,
            0
        )
    ),
    SeparationPower_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "SeparationPower",
            this->coeffDict_,
            1
        )
    ),
    SeparationXPower_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "SeparationXPower",
            this->coeffDict_,
            1
        )
    ),
    activationIIFac_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "activationIIFac",
            this->coeffDict_,
            0
        )
    ),
    A0I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0I",
            this->coeffDict_,
            0
        )
    ),
    A1I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A1I",
            this->coeffDict_,
            0
        )
    ),
    A2I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A2I",
            this->coeffDict_,
            0
        )
    ),
    T3activation_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "T3activation",
            this->coeffDict_,
            0.0
        )
    ),
    B0I_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "B0I",
            this->coeffDict_,
            0.0
        )
    ),
    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            this->coeffDict_,
            false
        )
    ),
    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kInf",
            this->coeffDict_,
            k_.dimensions(),
            0
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaInf",
            this->coeffDict_,
            omega_.dimensions(),
            0
        )
    ),
    bijDelta_
    (
        IOobject
        (
            IOobject::groupName("bijDelta", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*dev(symm(fvc::grad(U)))/(omega_ + this->omegaMin_)
    ),
    Rij_
    (
        IOobject
        (
            IOobject::groupName("Rij", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0.0*(((2.0/3.0)*I)*k_ - this->nut_*twoSymm(fvc::grad(this->U_)))
    ),
    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha2", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*k_/k_
    ),
    alpha3_
    (
        IOobject
        (
            IOobject::groupName("alpha3", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*k_/k_
    ),
    separationFactor_
    (
        IOobject
        (
            IOobject::groupName("separationFactor", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        0*k_/k_
    ),
    SecondaryFlow_relaxation_
    (   
        dimensioned<scalar>::getOrAddToDict
        (
            "SecondaryFlow_relaxation",
            this->coeffDict_,
            0.5
        )
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    setDecayControl(this->coeffDict_);

    Rij_ = ((2.0/3.0)*I)*k_ - this->nut_*twoSymm(fvc::grad(this->U_)) + 2*k_*bijDelta_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
void kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::setDecayControl
(
    const dictionary& dict
)
{
    decayControl_.readIfPresent("decayControl", dict);

    if (decayControl_)
    {
        kInf_.read(dict);
        omegaInf_.read(dict);

        Info<< "    Employing decay control with kInf:" << kInf_
            << " and omegaInf:" << omegaInf_ << endl;
    }
    else
    {
        kInf_.value() = 0;
        omegaInf_.value() = 0;
    }
}


template<class BasicEddyViscosityModel>
bool kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        setDecayControl(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicEddyViscosityModel>
void kOmegaSSTWindFarmBase<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    BasicEddyViscosityModel::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));


// New calculations

    volTensorField gradU = fvc::grad(U);
    volSymmTensorField Sij = symm(gradU);
    volTensorField Oij = -0.5*(gradU - gradU.T());
    volScalarField S = sqrt(2*magSqr(symm(fvc::grad(U))));

    volScalarField tScale = 1./max( S/a1_ + this->omegaMin_,omega_ + this->omegaMin_);
    volScalarField tScale2 = tScale*tScale;
    volScalarField tScale3 = tScale*tScale2;
    volScalarField tScale4 = tScale*tScale3;

    volScalarField i1_ = tScale2 * tr(Sij & Sij);
    volScalarField i2_ = tScale2 * tr(Oij & Oij);
    volScalarField i3_ = tScale3 * tr((Sij & Sij) & Sij);
    volScalarField i4_ = tScale3 * tr((Oij & Oij) & Sij);
    volScalarField i5_ = tScale4 * tr((Oij & Oij) & (Sij & Sij));

    volSymmTensorField T2_ = tScale2 * symm((Sij & Oij) - (Oij & Sij));
    volSymmTensorField T3_ = tScale2 * symm(Sij & Sij) - scalar(1.0/3.0)*I*i1_;


    forAll(Sij, CellI)
    {
        i1_[CellI] = tScale2[CellI] * tr(Sij[CellI] & Sij[CellI]);
        i2_[CellI] = tScale2[CellI] * tr(Oij[CellI] & Oij[CellI]);
        i3_[CellI] = tScale3[CellI] * tr((Sij[CellI] & Sij[CellI]) & Sij[CellI]);
        i4_[CellI] = tScale3[CellI] * tr((Oij[CellI] & Oij[CellI]) & Sij[CellI]);
        i5_[CellI] = tScale4[CellI] * tr((Oij[CellI] & Oij[CellI]) & (Sij[CellI] & Sij[CellI]));
        T2_[CellI] = tScale2[CellI] * symm((Sij[CellI] & Oij[CellI]) - (Oij[CellI] & Sij[CellI]));
        T3_[CellI] = tScale2[CellI] * symm(Sij[CellI] & Sij[CellI]) - scalar(1.0/3.0)*I*i1_[CellI];
    }


    

    if (SecondaryFlow_.value() == 1.0)
    {
        alpha2_ = -1.0*C0I_*k_/k_;
    }


    if (T3activation_.value() == 1.0)
    {
        alpha3_ = B0I_;
    }

    bijDelta_ = bijDelta_ + ((nut*omega_/(k_ + this->kMin_))*(alpha2_*T2_ + alpha3_*T3_) - bijDelta_)*SecondaryFlow_relaxation_;

    volSymmTensorField dAij = 2*k_*bijDelta_;
    volSymmTensorField P(-twoSymm(dAij & gradU));
    volScalarField Pk_bijDelta_ = 0.5*tr(P);

    dimensionedScalar nutMin("nutMin", dimensionSet(0, 2, -1, 0, 0, 0 ,0), 1e-9);

    volScalarField alpha1_separation = i1_ * 0.0;

    // Separation Factor

    if (Separation_.value() == 1.0)
    {
        alpha1_separation = -1.0*(A0I_*k_/k_);
        Info << "if clause is activated" << endl;
    }

    separationFactor_ = pow(
        max
        (
            min
            (
                (scalar(1) - pow(nut*omega_/k_, SeparationXPower_.value())),
                scalar(1)
            ),
            scalar(0)
        ),SeparationPower_.value())*(alpha1_separation); 



//  Continue with kOSST


    volScalarField::Internal GbyNu0
    (
        this->type() + ":GbyNu",
        // ((tgradU() && dev(twoSymm(tgradU()))) )
        ( (tgradU() && dev(twoSymm(tgradU()))) + Pk_bijDelta_/(nut+ nutMin) )
    );

    volScalarField::Internal G(this->GName(), nut*GbyNu0);
    // const fvMesh& mesh1 = this->mesh_;
    // volSymmTensorField& bijDelta = mesh1.lookupObjectRef<volSymmTensorField>("bijDelta");

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();


    

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());

    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        GbyNu0 = GbyNu(GbyNu0, F23(), S2());

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            // alpha()*rho()*gamma*GbyNu0
            alpha()*rho()*gamma*(GbyNu0+GbyNu0*separationFactor_)  // Cw1 F1 omega/k Pk
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          - fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + alpha()*rho()*beta*sqr(omegaInf_)
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(F1), k_)
     ==
        alpha()*rho()*Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), k_)
      + alpha()*rho()*betaStar_*omegaInf_*kInf_
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    tgradU.clear();

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(S2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
