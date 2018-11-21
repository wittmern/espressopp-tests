#include "LJ_include.hpp"
//#include <fenv.h>

int main(int argc, char *argv[]) {
  initMPIEnv(argc, argv);
//  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // integration steps, cutoff, skin and thermostat flag
  espressopp::longint steps = 5000;
  espressopp::real timestep = 0.005;
  espressopp::longint intervals = 500;

  espressopp::real rc = 2.31; // CG cutoff, Morse
  espressopp::real rca = 1.122462048309373; // AT cutoff (2^(1/6)), WCA
  espressopp::real skin = 0.4;

  espressopp::real gamma = 0.5;
  espressopp::real temp = 1.0;

  espressopp::real ex_size = 12.5;
  espressopp::real hy_size = 5.0;
  
  boost::shared_ptr<espressopp::System> system(new espressopp::System());
  system->comm = mpiWorld;

  // tabulated morse potential used for CG interactions
  std::string tabMorse = "pot-morse.txt";
  boost::shared_ptr<espressopp::interaction::Morse> 
    potMorseTab(new espressopp::interaction::Morse(0.105,2.4,rc,rc));  
  
  espressopp::real low = 0.005;
  espressopp::real high = 4.5;
  espressopp::real body = 2;
  int N = 512;
  if(system->comm->rank() == 0) write_tab_file(potMorseTab, tabMorse, N, low, high, body);
  system->comm->barrier();
  
  std::string input_file = "adress.espressopp";
  espressopp::real Lx, Ly, Lz;
  std::vector<espressopp::real> x, y, z;
  std::vector<int> type, q;
  std::vector<espressopp::real> vx, vy, vz, fx, fy, fz;
  std::vector<std::pair<int,int>> bondpairs;
  
  read_espresso_data(input_file,
                      Lx, Ly, Lz,
                      x, y, z,
                      type, q,
                      vx, vy, vz,
                      fx, fy, fz,
                      bondpairs);

  espressopp::longint num_particlesCG = 5001;
  espressopp::longint num_particles = (espressopp::longint)x.size();

  std::cout << "Setting up simulation ... \n";

  espressopp::real density = num_particles / (Lx * Ly * Lz);
  espressopp::Real3D size(Lx, Ly, Lz);

  boost::shared_ptr<espressopp::esutil::RNG> rng (new espressopp::esutil::RNG());
  boost::shared_ptr<espressopp::bc::BC> bc(new espressopp::bc::OrthorhombicBC(system->rng, size));
  system->rng = rng;
  system->bc = bc;
  system->setSkin(skin);

  espressopp::Int3D nodeGrid(2, 2, 1);
  espressopp::Int3D cellGrid(6, 6, 12);

  boost::shared_ptr<espressopp::storage::Storage> storage (new espressopp::storage::DomainDecompositionAdress(system, nodeGrid, cellGrid));
  system->storage = storage;

  boost::shared_ptr<espressopp::FixedTupleListAdress> ftpl(new espressopp::FixedTupleListAdress(system->storage));
  boost::shared_ptr<espressopp::FixedPairListAdress> 
    fpl(new espressopp::FixedPairListAdress(system->storage, ftpl));
  
  // prepare AT particles  
  for(espressopp::longint pidCG=0;pidCG<num_particlesCG;++pidCG) {
    espressopp::Real3D cmp(0.0);
    espressopp::Real3D cmv(0.0);
    std::vector<espressopp::Particle*> tupleAT;
    
      for(espressopp::longint j=0;j<4;++j) {
        espressopp::longint pidAT = pidCG*4+j;
        espressopp::Real3D pos(x[pidAT], y[pidAT], z[pidAT]);
        espressopp::Real3D vel(vx[pidAT], vy[pidAT], vz[pidAT]);
        espressopp::Real3D force(0.0);
        for (espressopp::longint i=0;i<3;++i) {
          cmp[i] += pos[i] / 4.0;
          cmv[i] += vel[i] / 4.0;
        }
      }
      espressopp::Particle* particleCG = system->storage->addParticle(pidCG+num_particles,cmp);
      if(particleCG != 0) {
        particleCG->setV(cmv);
        particleCG->setMass(4.0);
        particleCG->setF(espressopp::Real3D(0.0));
        particleCG->setType(0);
      }
    
    ftpl->add(pidCG+num_particles);
    
    for(espressopp::longint pidAT=0;pidAT<4;++pidAT) {
      espressopp::longint pid = pidCG*4+pidAT;
      espressopp::Real3D pos(x[pidAT], y[pidAT], z[pidAT]);
      espressopp::Real3D vel(vx[pidAT], vy[pidAT], vz[pidAT]);
      espressopp::Real3D force(fx[pidAT], fy[pidAT], fz[pidAT]);
      size_t type = 1;
      espressopp::real mass = 1.0;
      espressopp::Particle* particleAT = system->storage->addAdrATParticle(pid, pos, cmp);
      if(particleAT != 0) {
        tupleAT.push_back(particleAT);
        particleAT->setV(vel);
        particleAT->setF(force);
        particleAT->setType(type);
        particleAT->setMass(mass);
      }
      
      ftpl->add(pid);
    }
    ftpl->addTs(); // add tuples per CG particle
  }

  system->storage->setFixedTuplesAdress(ftpl);
  
  
  // add bonds  
  for(auto it=bondpairs.begin();it!=bondpairs.end();++it) {
    fpl->add(it->first,it->second);
  } 

  std::cout << "Added tuples and bonds. Now decomposing...\n";
  system->storage->decompose();

  bool rebuildVL = true;
  boost::shared_ptr<espressopp::VerletListAdress> 
    verletList(new espressopp::VerletListAdress(system, rc+skin, rc+skin,rebuildVL, ex_size, hy_size));
  verletList->setAdrCenter(18.42225, 18.42225, 18.42225);

  // non-bonded potentials
  // LJ Capped WCA between AT and tabulated Morse between CG particle
  boost::shared_ptr<espressopp::interaction::VerletListAdressInteractionTemplate<espressopp::interaction::LennardJonesCapped, espressopp::interaction::Tabulated>>
    interNB (new espressopp::interaction::VerletListAdressInteractionTemplate<espressopp::interaction::LennardJonesCapped, espressopp::interaction::Tabulated>(verletList, ftpl));

  // AT
  espressopp::real epsilon = 1.0;
  espressopp::real sigma = 1.0;
  boost::shared_ptr<espressopp::interaction::LennardJonesCapped>
    potWCA (new espressopp::interaction::LennardJonesCapped(epsilon, epsilon, rca, 0.27));
  potWCA->setShift(1.0);
 
  // CG 
  int itype = 2;
  boost::shared_ptr<espressopp::interaction::Tabulated>
    potMorse (new espressopp::interaction::Tabulated(itype, tabMorse.c_str(), rc));
  
  interNB->setPotentialAT(1,1, *potWCA);
  interNB->setPotentialCG(0,0, *potMorse);

  system->addInteraction(interNB);

  //bonded potentials
  // FENE and LJ potential between AT particles
  espressopp::real K = 30.0;
  espressopp::real r0 = 0.0;
  espressopp::real rMax = 1.5;

  boost::shared_ptr<espressopp::interaction::FENE> 
    potFENE (new espressopp::interaction::FENE(K, r0, rMax, espressopp::infinity));

  boost::shared_ptr<espressopp::interaction::LennardJones> 
    potLJ (new espressopp::interaction::LennardJones(epsilon, sigma, rca));
  potLJ->setShift(1.0);

  boost::shared_ptr<espressopp::interaction::FixedPairListInteractionTemplate<espressopp::interaction::FENE>> 
    interFENE (new espressopp::interaction::FixedPairListInteractionTemplate<espressopp::interaction::FENE>(system, fpl, potFENE));

  boost::shared_ptr<espressopp::interaction::FixedPairListInteractionTemplate<espressopp::interaction::LennardJones>>
    interLJ (new espressopp::interaction::FixedPairListInteractionTemplate<espressopp::interaction::LennardJones>(system, fpl, potLJ));

  system->addInteraction(interFENE);
  system->addInteraction(interLJ);

  // VV integrator
  boost::shared_ptr<espressopp::integrator::VelocityVerlet>
    integrator (new espressopp::integrator::VelocityVerlet(system));
  integrator->setTimeStep(timestep);

  // add AdResS extension
  boost::shared_ptr<espressopp::integrator::Adress> 
    adress (new espressopp::integrator::Adress (system, verletList, ftpl));
  integrator->addExtension(adress);

  // add Langevin thermostat extension
  boost::shared_ptr<espressopp::integrator::LangevinThermostat>
    langevin (new espressopp::integrator::LangevinThermostat(system));
  langevin->setGamma(gamma);
  langevin->setTemperature(temp);
  langevin->setAdress(true);
  integrator->addExtension(langevin);

  if(system->comm->rank() == 0) {
  std::cout << "\n" 
            << "number of AT particles = " << num_particles << "\n"
            << "number of CG particles = " << num_particlesCG << "\n"
            << "density = " << std::setprecision(4) << density << "\n"
            << "rc = " << rc << "\n"
            << "dt = " << integrator->getTimeStep() << "\n"
            << "skin = " << system->getSkin() << "\n"
            << "steps = " << steps << "\n"
            << "NodeGrid = " << nodeGrid << "\n"
            << "CellGrid = " << cellGrid << "\n"
            << "\n";
  }
  
  // Initizalize storage dump
  std::string dumpFile = "testfile.dat";
  bool unfolded = true;
  espressopp::real lengthFactor = 1;
  std::string lengthUnit = "LJ";
  bool storePIDs = true;
  bool storeVelocities = true;
  bool append = false;

  system->comm->barrier();
  
  espressopp::io::DumpXYZ xyzDumper(system,integrator, "data.xyz", unfolded, lengthFactor, lengthUnit, storePIDs, storeVelocities, append);
  xyzDumper.dump();

  // analysis
  boost::shared_ptr<espressopp::analysis::Temperature>
    temperature (new espressopp::analysis::Temperature(system));
  boost::shared_ptr<espressopp::analysis::Pressure>
    pressure (new espressopp::analysis::Pressure(system));
  boost::shared_ptr<espressopp::analysis::PressureTensor>
    pressureTensor (new espressopp::analysis::PressureTensor(system));
  
  espressopp::real T = temperature->compute_real();
  espressopp::real P = pressure->compute();
  espressopp::Tensor Pij = pressureTensor->computeRaw();
  espressopp::real Ek = T * 0.5 * (3 * num_particles);
  espressopp::real Ep = interNB->computeEnergy();
  espressopp::real Eb = interFENE->computeEnergy();
  if(system->comm->rank() == 0) {
    std::cout << " step     T          P        Pxy       etotal     epotential      ebonded     ekinetic\n";
    std::cout << std::setw(5) << 0 
              << std::setw(8) << std::setprecision(4) << T 
              << std::setw(10) << std::setprecision(5) << P
              << std::setw(8) << std::setprecision(5) << Pij[3]
              << std::setw(12) << std::setprecision(3) << Ek + Eb  + Ep
              << std::setw(12) << std::setprecision(3) << Ep 
              << std::setw(12) << std::setprecision(3) << Eb 
              << std::setw(12) << std::setprecision(3) << Ek 
              << "\n";
  }
  system->comm->barrier();
  time_t timer;
  time(&timer);
  espressopp::real nsteps = steps / intervals;
  for (espressopp::longint s=1;s<intervals+1;++s) {
    integrator->run(nsteps);
    xyzDumper.dump();
    espressopp::longint step = nsteps * s;
    T = temperature->compute_real();
    P = pressure->compute();
    Pij = pressureTensor->computeRaw();
    Ek = 0.5 * T * (3 * num_particles);
    Eb = interFENE->computeEnergy();
    Ep = interNB->computeEnergy();
    if(system->comm->rank() == 0) {
    std::cout << std::setw(5) << step 
              << std::setw(8) << std::setprecision(4) << T 
              << std::setw(10) << std::setprecision(5) << P 
              << std::setw(8) << std::setprecision(5) << Pij[3] 
              << std::setw(12) << std::setprecision(3) << Ek + Eb  + Ep
              << std::setw(12) << std::setprecision(3) << Ep 
              << std::setw(12) << std::setprecision(3) << Eb 
              << std::setw(12) << std::setprecision(3) << Ek 
              << "\n";
    }
    system->comm->barrier();
    system->storage->decompose();
  }
  double seconds = difftime(timer,time(NULL));

  if(system->comm->rank() == 0) {
  // simulation information
  std::cout << "Neighbor list builds = " << verletList->getBuilds() << "\n"
            << "Integration steps = " << integrator->getStep() << "\n"
            << "CPU time = " << std::setprecision(1) << seconds << "\n";
  }
  
  finalizeMPIEnv();


	return 0;
}
