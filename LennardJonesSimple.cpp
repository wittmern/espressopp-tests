#include <cmath>
#include <memory>

#include "System.hpp"
#include "Real3D.hpp"
#include "Int3D.hpp"
#include "include/mpi.hpp"
#include "VerletList.hpp"
#include "interaction/LennardJones.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "storage/DomainDecomposition.hpp"
#include "storage/CellGrid.hpp"
#include "esutil/RNG.hpp"
#include "interaction/VerletListInteractionTemplate.hpp"
#include "include/esconfig.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "storage/NodeGrid.hpp"
#include "main/espressopp_common.hpp"
#include "boost/make_shared.hpp"
#include "io/DumpXYZ.hpp"


int main(int argc, char *argv[]) {
  initMPIEnv(argc, argv);

  
  boost::shared_ptr<class espressopp::System> system (new espressopp::System());
  
  boost::shared_ptr<espressopp::esutil::RNG> rng (new espressopp::esutil::RNG());
  const espressopp::Real3D box(30, 30, 30);
  
  boost::shared_ptr<espressopp::bc::BC> bc (new espressopp::bc::OrthorhombicBC(rng, box));

  system->bc   = bc;
  system->rng  = rng;
  system->comm = mpiWorld;

  espressopp::real maxcutoff = pow(2,3);
  espressopp::real skin = 0.4;

  system->setSkin(skin);

  // Number of nodes
  espressopp::Int3D nodegrid_size(2,1,1);
  // Number of Cells
  espressopp::Int3D cellgrid_size(12,12,12);

  //initialize Storage
  boost::shared_ptr<espressopp::storage::Storage> ddStorage (new espressopp::storage::DomainDecomposition(system->getShared(),nodegrid_size,cellgrid_size) );
  system->storage = ddStorage;

  //initialize integrator
  boost::shared_ptr<espressopp::integrator::VelocityVerlet> 
    integrator (new espressopp::integrator::VelocityVerlet(system));
  
  espressopp::real dt = 0.005;
  integrator->setTimeStep(dt);
  
  //add some particle to the system
  espressopp::longint i=0;
  for(espressopp::real x=0;x<26;++x) {
    for(espressopp::real y=0;y<26;++y) {
      for(espressopp::real z=0;z<26;++z) {
        espressopp::Real3D pos(
                            x*30.0/26.0,
                            y*30.0/26.0,
                            z*30.0/26.0);
        ddStorage->addParticle(i,pos);
        i++;
      }}}

  boost::shared_ptr<espressopp::VerletList> verletList( new espressopp::VerletList(system, maxcutoff, true));

  // initizalize potential
  boost::shared_ptr<espressopp::interaction::LennardJones> 
    LJPot(new espressopp::interaction::LennardJones(1.0,1.0,maxcutoff, 1));

  // initialize interaction using potential and integrator
  boost::shared_ptr<espressopp::interaction::VerletListInteractionTemplate<espressopp::interaction::LennardJones>>  
    nonBondedInteraction (new espressopp::interaction::VerletListInteractionTemplate<espressopp::interaction::LennardJones>(verletList));
  
  nonBondedInteraction->setPotential(0,0,*LJPot);
  system->addInteraction(nonBondedInteraction);

  
  // Initizalize storage dump
  std::string dumpFile = "testfile.dat";
  bool unfolded = false;
  espressopp::real lengthFactor = 1;
  std::string lengthUnit = "LJ";
  bool storePIDs = true;
  bool storeVelocities = true;
  bool append = false;
  espressopp::io::DumpXYZ dumper(system, integrator, dumpFile, unfolded, lengthFactor, lengthUnit, storePIDs, storeVelocities, append);
  dumper.dump();

  std::cout << ddStorage->getInt3DCellGrid() << std::endl;

  // run integration and dump particle every ten steps
  for (int k=0;k<10;++k) {
    integrator->run(10);
    dumper.dump();
    std::cout << "Finish Step " << k*10 << "..." << std::endl;
  }

  finalizeMPIEnv();


	return 0;
}


