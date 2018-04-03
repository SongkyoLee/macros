int run_macro( 
	      std::string infile = "/sphenix/user/dvp/sims/IHCal/jet/output_00000.root",
	      std::string outfile = "test.root"
	       )
{
  
  gSystem->Load("libg4dst.so");
  gSystem->Load("libfun4all.so");
  gSystem->Load("libphfield_io.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libphhepmc.so");
  gSystem->Load("libg4testbench.so");
  gSystem->Load("libg4hough.so");
  gSystem->Load("libcemc.so");
  gSystem->Load("libg4eval.so");
  gSystem->Load("libcalotrigger.so");

  gSystem->Load("libJetTreeMaker.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(20);

  recoConsts *rc = recoConsts::instance();

  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  hitsin->fileopen( infile );
  se->registerInputManager(hitsin);

  JetTreeMaker *tm = new JetTreeMaker( outfile );
  se->registerSubsystem( tm );

  se->run();

  se->End();
  std::cout << "All done" << std::endl;
  delete se;

  gSystem->Exit(0);
}
