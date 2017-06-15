classdef Test_compPennation < matlab.unittest.TestCase
    
    methods(Test)
	
%% unit tests for compPennation
       
        function test_compPennation_tetrahedron(testCase)
			load('+tests/compPennationTestData.mat')
            [testSTL.vertices,testSTL.faces, testSTL.normals] = stlread('+tests/test.STL');
            %note: order of structure fields depends on version of stlread
			
            [actSolution1, actSolution2] = compPennation(testSTL,origin,insertion,mean_vector);
            expSolution1 = 44.311707716769803;
            expSolution2 = 44.311707716769803;
            verifyEqual(testCase,actSolution1,expSolution1);
            verifyEqual(testCase,actSolution2,expSolution2);
        end
                       
    end
    
end

