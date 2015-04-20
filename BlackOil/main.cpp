#include "ForwardSimulation.hpp"
#include "Adjoint.hpp"
#include "CSR_Example.hpp"





int main ()
{
   Simulator sim;
   sim.Run();


   ADJOINT adjoint;

   std::vector< std::vector<double> > gi_dm_i;
   adjoint.Gi_Dm_i( sim, sim.Get_RecordWellData()[1], gi_dm_i );




   

   //int TSIndex = sim.GetTimeStepIndex_RecordWellData().back();
   // double TS = sim.GetTimeStep()[ TSIndex ];
   // double T = sim.GetTotalTime()[ TSIndex ];

   // std::cout << "Time Step = " << TS << std::endl;
   // std::cout << "Time      = " << T  << std::endl;

   // const std::vector<double>& record = sim.GetRecord()[ TS ];
   // std::cout << "Recorded Data"<< std::endl;
   // for( std::size_t i = 0; i < record.size(); ++i )
   // {
   //    std::cout << i << "\t" << record[i] << std::endl;
   // }
   
   // ADJOINT adjoint;
   // CSR<> J, Jm;
   // std::vector<double> gi_dy_col;
   // adjoint.RecoverJacobian( record, TS, J );
   // adjoint.Gi_Dy_Col_Producer( sim.GetProblem(),
   // 			       sim.GetProblem().GetProducer()[0],
   // 			       record,
   // 			       gi_dy_col );
   // std::vector<double> lamda_col;
   // adjoint.Lamda_Col( J, gi_dy_col, lamda_col );
			       
   // adjoint.RecoverJacobian_Dm( record, TS, Jm );

   // std::vector<double> result = lamda_col * Jm;
   // std::cout << "result"<< std::endl;
   // for( std::size_t i = 0; i < result.size(); ++i )
   // {
   //    std::cout << i << "\t" << result[i] << std::endl;
   // }
   
   

   

   return -1;
}
