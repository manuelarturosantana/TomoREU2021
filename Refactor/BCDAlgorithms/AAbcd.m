function [x,BCDinfo] = AAbcd(x_0,b,options,probinfo)
    %Intialization for Anderson Acceloration. Should come from something passed
    %in probably
    x_curr = x0;
    BCDinfo.x = x_curr;
    G = [ ];
    num_stored_residuals = 0;
    %Number of previous solutions to use with Anderson Acceloration


    %This enters the AABCD optimization loop.
    for k = 1:optIter
        [g_curr, iterInfo] = fpBCD(BCDinfo);
        f_curr = g_curr - x_curr;
        if k > 1
            %Form G on the second iteration.
            delta_f = f_curr - f_old;
            delta_g = g_curr - g_old;

            %Update the columns on G
            if num_stored_residuals < max_stored_residuals
                G = [G delta_g];
            else
                G = [G(:,2:end) delta_g];
            end
            %We update this even though the second step doesn't actually
            %increase the number of columns to G. This is a work saving step
            %because we have to redo the QR factorization of F, but will need G
            %in this size later.
            num_stored_residuals = num_stored_residuals + 1;
        end
        %Set f and g old for calculating delta f the next iteration.
        f_old = f_curr; g_old = g_curr;
        %Sets the x value on the first iteration
        if num_stored_residuals == 0
            x_curr = g_curr;
        else
            %one vector QR factorization
            if num_stored_residuals == 1 
                Q(:,1) = delta_f / norm(delta_f);
                R(1,1) = norm(delta_f);
            else
                %delete a QR column if too many.
                if num_stored_residuals > max_stored_residuals 
                    [Q,R] = qrdelete(Q,R,1);
                    num_stored_residuals = num_stored_residuals - 1;
                    %This deals with a qrdelete usage explained below.
                    if size(R,1) ~= size(R,2)
                        Q = Q(:,1:num_stored_residuals -1);
                        R = R(1:num_stored_residuals - 1,:);
                    end
                    %Explination: If Q is not square then the matlab function
                    %QR delete removes one column of Q, and and column and row
                    %of R. If Q is square then the column demenion of Q is not
                    %reduced, and R only has a column removed. This behavior is
                    %to account for thick QR decompoitions, but since we are 
                    %using a thin one we must correct if this happens.
                end
                %One pass of modified Gram-Schmidt to update the QR
                %factorization. Recall this uses the fact that since F = QR
                %then the column of F is a linear combination of columns of Q
                %with coefficents coming from the last column of R.
                for i = 1:num_stored_residuals -1 
                    R(i,num_stored_residuals) = Q(:,i)' * delta_f;
                    delta_f = delta_f - R(i,num_stored_residuals) * Q(:,i); 
                end
                %Completing the Graham-Scmidt Iteration
                R(num_stored_residuals,num_stored_residuals) = norm(delta_f);
                Q = [Q,delta_f / norm(delta_f)]; 
            end
            %Here we delete more columns in the QR factorization to deal with
            %poor conditioning of the R matrix that may occur.
            while (cond(R)) > dropTol && num_stored_residuals > 1
               [Q,R] = qrdelete(Q,R,1);
               num_stored_residuals = num_stored_residuals - 1;
               %If we change the size of the QR, we must also change the size
               % of G.
               G = G(2:end);
               %See above explination for this step
               if size(R,1) ~= size(R,2)
                        Q = Q(:,1:num_stored_residuals -1);
                        R = R(1:num_stored_residuals - 1,:);
               end
            end
           gamma = R \ (Q' * f_curr);
           %updating the x info based on the anderson acceleration.
           x_curr = g_curr - G * gamma;
        end

    %Update the current x guess in the BCD info.
    BCDinfo.x = x_curr;
    BCDinfo.RParams = iterInfo.RParams;
    BCDinfo.angleParams = iterInfo.angleParams;
    %Store the errors for plotting at each iteration.
    xs = [xs, x_curr];
    xErrors = [xErrors,norm(x_curr - xtrue) / norm(xtrue)];
    pErrors = [pErrors,norm(paramTrue - iterInfo.p_0)/norm(paramTrue)];
    RErrors = [RErrors,norm(RPert - iterInfo.RParams) / norm(RPert)];
    angErrors = [angErrors,norm(iterInfo.angleParams - angle_pert)/norm(angle_pert)];
    end    

end