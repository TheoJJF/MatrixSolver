classdef DenseValidation < matlab.unittest.TestCase

    methods (Test)

        function case1(testCase)
            A = [2,1,1;0,3,2;0,0,4];
            b = [7;12;12];

            actual = SolveGEPP(A,b);
            expected = A\b;
            testCase.verifyEqual(actual, expected); 
        end

        function case2(testCase)
            A = [2,4,6;1,5,-3;3,1,8];
            b = [12;13;12];

            actual = SolveGEPP(A,b);
            expected = A\b;
            testCase.verifyEqual(actual, expected); 
        end

        function case3(testCase)
            A = [1,1,0;2,3,1;4,1,1];
            b = [1;1;2];

            actual = SolveGEPP(A,b);
            expected = A\b;
            testCase.verifyEqual(actual, expected); 
        end

        function case4(testCase)
            A = [7,2,-3;3,1,8;10,3,5];
            b = ones(3,1);

            testCase.verifyError(@()SolveGEPP(A,b),"GEPP:SingularMatrix"); 
        end
    end

end