#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union
from collections import namedtuple
import sys
import os

from acts.examples import Sequencer, GenericDetector, RootParticleReader

import argparse
import pathlib
import acts
import acts.examples

from acts import UnitConstants as u

from acts.examples.simulation import (
    addParticleGun,
    addParticleReader,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    ParticleSmearingSigmas,
    addCKFTracks,
    TrackSelectorConfig,
    # CKFPerformanceConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    CkfConfig,
    TruthEstimatedSeedingAlgorithmConfigArg,
)
from acts.examples import TGeoDetector




def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(
        description="Command line arguments for CKF")
    parser.add_argument(
        "-i",
        "--indir",
        dest="indir",
        help="Directory with input root files",
        default="./",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="outdir",
        help="Output directory for new ntuples",
        default=None,
    )
    parser.add_argument(
        "-n", "--nEvents", dest="nEvts", help="Number of events to run over", default=1, type=int
    )

    parser.add_argument(
        '-g',
        '--geant4',
        help='Use Geant4',
        action='store_true'
    )
    
    parser.add_argument(
        '-d',
        '--dead',
        help='Apply dead zones',
        action='store_true'
    )
    parser.add_argument(
        '-f',
        '--fast',
        help='Apply fast sim selections',
        action='store_true'
    )

    parser.add_argument(
        '-p',
        '--gun',
        help='Use particle gun',
        type=int,
        default=0,
    )

    parser.add_argument(
        '-s',
        '--truthSeeding',
        help='Use particle truth to perform the seeding',
        action='store_true'
    )
    parser.add_argument(
        '-m',
        '--truthSmeared',
        help='Use particle truth smeared to perform the seeding',
        action='store_true'
    )
    parser.add_argument(
        '-v',
        '--truthVertexing',
        help='Use particle truth to perform the vertexing',
        action='store_true'
    )
    parser.add_argument(
        "--sf_maxSeedsPerSpM",
        dest="sf_maxSeedsPerSpM",
        help="Number of compatible seeds considered for middle seed",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--sf_cotThetaMax",
        dest="sf_cotThetaMax",
        help="cot of maximum theta angle",
        type=float,
        default=7.40627,
    )
    parser.add_argument(
        "--sf_sigmaScattering",
        dest="sf_sigmaScattering",
        help="How many sigmas of scattering to include in seeds",
        type=float,
        default=5,
    )
    parser.add_argument(
        "--sf_radLengthPerSeed",
        dest="sf_radLengthPerSeed",
        help="Average Radiation Length",
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--sf_impactMax",
        dest="sf_impactMax",
        help="max impact parameter in mm",
        type=float,
        default=3.0,
    )
    parser.add_argument(
        "--sf_maxPtScattering",
        dest="sf_maxPtScattering",
        help="maximum Pt for scattering cut in GeV",
        type=float,
        default=10.0,
    )
    parser.add_argument(
        "--sf_deltaRMin",
        dest="sf_deltaRMin",
        help="minimum value for deltaR separation in mm",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--sf_deltaRMax",
        dest="sf_deltaRMax",
        help="maximum value for deltaR separation in mm",
        type=float,
        default=60.0,
    )

    return parser


def runCKFTracks(
    detector,
    trackingGeometry,
    inputDir: Path,
    outputDir: Path,
    useGeant4=False,
    truthSeeding=False,
    truthSmeared=False,
    truthVertexing=False,
    applyDeadZones=False,
    applyFastSimSelections=False,
    particleGun=0,
    NumEvents=1,
    MaxSeedsPerSpM=3,
    CotThetaMax=5.809222379141632,
    SigmaScattering=7.293572849910582,
    RadLengthPerSeed=0.07827007716884793,
    ImpactMax=10.47494916604592,
    MaxPtScattering=10.0
    # DeltaRMin=1.0,
    # DeltaRMax=60.0,
):

    field = acts.ConstantBField(acts.Vector3(
        0.0, 1.5 * u.T, 0.0))  # ~dipole field

    rnd = acts.examples.RandomNumbers(seed=44)

    s = acts.examples.Sequencer(
        events=NumEvents,
        # events=1,
        numThreads=1,
        outputDir=str(outputDir),
        # to remove mechanism to catch FPE during algorithm execution (see Mattermost 18/7/23 Igor)
        trackFpes=False,
    )

    if particleGun >0:
        addParticleGun(
            s,
            MomentumConfig(0 * u.GeV, 100* u.GeV, transverse=False),
            EtaConfig(0,5.2, uniform=True),
            ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
            multiplicity=particleGun,
            vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(0.0 * u.mm, 0.0 * u.mm, 0.0 * u.mm, 0 * u.ns),
            mean=acts.Vector4(0, 0, 0, 0),
            ),

            rnd=rnd,
            outputDirRoot=outputDir,
        )
    else:
        addParticleReader(
            s,
            inputDir=inputDir,
            outputDirRoot=outputDir,
        )
    

    if useGeant4:
        addGeant4(
                    s,
                    detector,
                    trackingGeometry,
                    field,
                    rnd,
                    #g4DetectorConstructionFactory: Optional[Any] = None,
                    #volumeMappings: List[str] = [],
                    #materialMappings: List[str] = [],
                    preSelectParticles=ParticleSelectorConfig(
                        eta=(None, None),
                        pt=(1000 * u.MeV, None),
                        removeNeutral=True,
                    ),
                    #recordHitsOfSecondaries=True,
                    #keepParticlesWithoutHits=True,
                    outputDirRoot=outputDir,
                    killVolume=trackingGeometry.worldVolume,
                    killAfterTime=25 * u.ns,
                )
    else:
        addFatras(
            s,
            trackingGeometry,
            field,
            rnd,
            preSelectParticles=ParticleSelectorConfig(
                eta=(None, None),
                pt=(100 * u.MeV, None),
                removeNeutral=True,
            ),
            outputDirRoot=outputDir,
            # logLevel=acts.logging.DEBUG,
        )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile="geomVTNA60p/digismear.json",
        outputDirRoot=outputDir,
        efficiency=0.99 if applyDeadZones or applyFastSimSelections else 1,
        applyDeadAreas=applyDeadZones,
        applyFastSimSelections=applyFastSimSelections,
        rnd=rnd,
    )

    # rotation of B field, so that B = Bz, needed for seeding
    field2 = acts.ConstantBField(acts.Vector3(
        0.0, 0.0, -1.5 * u.T))  # ~dipole field

    addSeeding(
        s,
        trackingGeometry,
        field2,
        TruthSeedRanges(pt=(100 * u.MeV, None),
                        eta=(None, None), nHits=(4, None)),
        SeedFinderConfigArg(
            maxSeedsPerSpM=MaxSeedsPerSpM,
            cotThetaMax=CotThetaMax,
            sigmaScattering=SigmaScattering,
            radLengthPerSeed=RadLengthPerSeed,
            maxPtScattering=MaxPtScattering,
            # min and max R between Middle and Top SP
            deltaRTopSP=(45 * u.mm, 100 * u.mm),
            # min and max R between Middle and Bottom SP
            deltaRBottomSP=(50 * u.mm, 100 * u.mm),
            impactMax=ImpactMax * u.mm,
            deltaZMax=50,  # was 5
            minPt=100 * u.MeV,
            interactionPointCut=True,
            useVariableMiddleSPRange=False,  # MODIFICATO 22/5/23
            # not useful if useVariableMiddleSPRange=False
            deltaRMiddleSPRange=(0 * u.mm, 0 * u.mm),
            collisionRegion=(-0.5 * u.mm, 0.5 * u.mm),  # 0.5
            # NOT USED??? seems to be used in Orthogonal seeding
            r=(0 * u.mm, 400 * u.mm),
            rMiddle=(149 * u.mm, (250**2+150**2/2.)**(0.5) * u.mm),
            # deltaR=(0.01 * u.mm, 2.5 * u.mm),      #NOT USED???
            # z=(20 * u.mm, 80 * u.mm),      #NOT USED??? seems to be used in Orthogonal seeding
            # zOutermostLayers={74,76}, #modified by me in SeedfinderConfig.hpp
            seedConfirmation=False,
            forwardSeedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=150 * u.mm,
                zMaxSeedConf=0 * u.mm,
                rMaxSeedConf=7 * u.mm,
                nTopForLargeR=1,
                nTopForSmallR=1,
            ),
            centralSeedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=7 * u.mm,
                zMaxSeedConf=-7 * u.mm,
                rMaxSeedConf=80 * u.mm,
                nTopForLargeR=1,
                nTopForSmallR=1,
            ),
        ),
        SeedFinderOptionsArg(
            beamPos=(0 * u.mm, 0 * u.mm),
            bFieldInZ=1.5 * u.T,
        ),  # why should I give the b field? to compute phi bins in SpacePointGrid.ipp
        SeedFilterConfigArg(  # not used, why?
            seedConfirmation=False,
            maxSeedsPerSpMConf=5,
            maxQualitySeedsPerSpMConf=5,
            compatSeedLimit=2,  # added 4/10/23 (value not tuned)
        ),
        SpacePointGridConfigArg(
            rMax=400 * u.mm,
            zBinEdges=[-160., -160., 160., 160.],
            # not used if bfieldZ is 0, otherwise it's used to compute number of phi bins
            impactMax=0.1 * u.mm,
            # not used if bfieldZ is 0, otherwise it's used to compute number of phi bins
            phiBinDeflectionCoverage=1,
        ),
        SeedingAlgorithmConfigArg(
            zBinNeighborsTop=[[0, 1], [-1, 1], [-1, 0]],
            zBinNeighborsBottom=[[0, 1], [-1, 1], [-1, 0]],
            numPhiNeighbors=1,
        ),
        TruthEstimatedSeedingAlgorithmConfigArg(deltaR=(0 * u.mm, 100000 * u.mm)),
        seedingAlgorithm=SeedingAlgorithm.TruthEstimated if truthSeeding else (SeedingAlgorithm.TruthSmeared if truthSmeared else SeedingAlgorithm.Default),
        geoSelectionConfigFile="geomVTNA60p/seed_config.json",
        outputDirRoot=outputDir,
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CkfConfig(
            chi2CutOff=30000,
            numMeasurementsCutOff=1,
            maxSteps=None,
        ),
        TrackSelectorConfig(
            pt=(100 * u.MeV, None),
            absEta=(None, None),
            nMeasurementsMin=4,
        ),
        outputDirRoot=outputDir,
    )
    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
        maximumSharedHits=1,
        nMeasurementsMin=4,
        maximumIterations=1000,
        ),
        outputDirRoot=outputDir,
    )

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Truth if truthVertexing else VertexFinder.Iterative,
        outputDirRoot=outputDir,
    )

    return s


if "__main__" == __name__:
    options = getArgumentParser().parse_args()

    # "events_full_sim"#"fullchain_input"
    inputDir = pathlib.Path.cwd() / "events_40GeV"
    suffix = ""
    if options.truthSeeding:
        suffix += "_truthEstimated"
    elif options.truthSmeared:
        suffix += "_particleSmearing"
    else:
        suffix += "_standardSeeding"
    if options.truthVertexing:
        suffix += "_truthVertexing"
    else:
        suffix += "_iterativeVertexing"
    if options.geant4:
        suffix += "_geant4"
    if options.gun > 0:
        suffix += "_gun"+str(options.gun)
    if options.dead:
        suffix += "_deadZones"
    current_dir = pathlib.Path.cwd()
    outputDir = str(current_dir / ("output" + suffix))

    matDeco = acts.IMaterialDecorator.fromFile(
        "geomVTNA60p/material-map_VTNA60p.json")
    jsonFile = "geomVTNA60p/tgeo-config_VTNA60p.json"
    tgeo_fileName = "geomVTNA60p/geom_VTNA60p.root"

    logLevel = acts.logging.INFO
    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)

    detector, trackingGeometry, decorators = TGeoDetector.create(jsonFile=str(jsonFile),
                                                                 fileName=str(
        tgeo_fileName),
        surfaceLogLevel=customLogLevel(),
        layerLogLevel=customLogLevel(),
        volumeLogLevel=customLogLevel(),
        mdecorator=matDeco,
    )
    print(outputDir if options.outdir == None else options.outdir)
    runCKFTracks(
        detector,
        trackingGeometry,
        inputDir=inputDir,
        outputDir=outputDir if options.outdir == None else options.outdir,
        useGeant4=options.geant4,
        particleGun=options.gun,
        applyDeadZones=options.dead,
        applyFastSimSelections=options.fast,
        truthSeeding=options.truthSeeding,
        truthSmeared=options.truthSmeared,
        truthVertexing=options.truthVertexing,
        NumEvents=options.nEvts,
        MaxSeedsPerSpM=options.sf_maxSeedsPerSpM,
        CotThetaMax=options.sf_cotThetaMax,
        SigmaScattering=options.sf_sigmaScattering,
        RadLengthPerSeed=options.sf_radLengthPerSeed,
        ImpactMax=options.sf_impactMax,
        # MaxPtScattering=options.sf_maxPtScattering,
        # DeltaRMin=options.sf_deltaRMin,
        # DeltaRMax=options.sf_deltaRMax,
    ).run()

