classdef TestHarness < matlab.unittest.TestCase
    
    methods(Test)
	
%% unit tests for YiCao_knnsearch
        function testYiCao_knnsearch_Basic(testCase)
            actSolution = YiCao_knnsearch(1,1,1);
            expSolution = 1;
            verifyEqual(testCase,actSolution,expSolution);
        end
        
        function testYiCao_knnsearch_KisTwoError(testCase)
            verifyError(testCase,@()YiCao_knnsearch(1,1,2),'MATLAB:badsubscript');
        end
        
        function testYiCao_knnsearch_KisTwo(testCase)
            actSolution = YiCao_knnsearch([1;1],[1;1],2);
            expSolution = [2 1;1 2];
            verifyEqual(testCase,actSolution,expSolution);
        end
        
        function testYiCao_knnsearch_multiOut(testCase)
            [actSolution1,actSolution2] = YiCao_knnsearch([1;1],[1;1],1);
            expSolution1 = [2;1];
            expSolution2 = [0;0];
            verifyEqual(testCase,actSolution1,expSolution1);
            verifyEqual(testCase,actSolution2,expSolution2);
        end
        
        function testYiCao_knnsearch_multiOut_Coords_BigRef(testCase)
            [actSolution1,actSolution2] = YiCao_knnsearch([1 1 1;2 2 2;3 3 3;4 4 4],[1 1 1; 3 3 3],1);
            expSolution1 = [1;3];
            expSolution2 = [0;0];
            verifyEqual(testCase,actSolution1,expSolution1);
            verifyEqual(testCase,actSolution2,expSolution2);
        end
		
%% unit tests for ipdm.m
		function test_ipdm_OneDataSet(testCase)
            actSolution = ipdm([1 1 1]);
            expSolution = 0;
            verifyEqual(testCase,actSolution,expSolution);
        end
		
		function test_ipdm_OneDataSet_multiRows(testCase)
            actSolution = ipdm([1 1 1;2 2 2; 3 3 3]);
            expSolution = [0 sqrt(3) sqrt(12);...
						sqrt(3) 0 sqrt(3); ...
						sqrt(12) sqrt(3) 0];
            verifyEqual(testCase,actSolution,expSolution);
        end
		
		function test_ipdm_TwoDataSets(testCase)
            actSolution = ipdm([1 1 1;2 2 2; 3 3 3],[1 1 1;2 2 2; 3 3 3]);
            expSolution = [0 sqrt(3) sqrt(12);...
						sqrt(3) 0 sqrt(3); ...
						sqrt(12) sqrt(3) 0];
            verifyEqual(testCase,actSolution,expSolution);
        end
		
		function test_ipdm_TwoDataSets_MismatchSize(testCase)
            actSolution = ipdm([1 1 1;2 2 2; 3 3 3],[1 1 1;2 2 2; 3 3 3;2 2 2; 1 1 1]);
            expSolution = [0 sqrt(3) sqrt(12) sqrt(3) 0;...
						sqrt(3) 0 sqrt(3) 0 sqrt(3); ...
						sqrt(12) sqrt(3) 0 sqrt(3) sqrt(12)];
            verifyEqual(testCase,actSolution,expSolution);
        end
                
    end
    
end

