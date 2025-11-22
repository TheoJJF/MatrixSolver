classdef SparseValidation < matlab.unittest.TestCase

    methods (Test)
        
        function case1(testCase)
            data = load(fullfile(pwd,"../data","SparseValidation.mat"));
            verification = load(fullfile(pwd,"../data","expected.mat"));

            A = data.A;
            b = data.b;
            expected = verification.expected;

            actual = SGEPP(A,b);
            testCase.verifyEqual(actual,expected);
        end

        

    end

end