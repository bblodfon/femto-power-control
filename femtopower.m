% NTUA Thesis code
% John Zobolas, May 2013

function femtopower
    format long;
    global results;
    results = cell(1,10);
    
    % the "standard" coordinates of the FAPs
    x = [200 400 400 700 600 600 300 850 800 200 600 800 500 200 100 400 700 900];
    y = [500 300 600 600 800 200 800 400 800 200 400 200 900 700 350 100 300 550];
    % Uncomment the below to get a random placement of the FAPs inside the Macrocell
    %x = randi([150 850],1,30);
    %y = randi([150 850],1,30);
    %length(x) = 18;
    xf = zeros(1,36);
    yf = zeros(1,36);
    for i=1:18
       while true
           random = 2*randi([-18 18],1,4);
            if all(random) % no zeros
                break;
            end
       end
       xf(2*i-1) = x(i)+random(1);% NRT users
       yf(2*i-1) = y(i)+random(2);
       xf(2*i) = x(i)+random(3);% RT users
       yf(2*i) = y(i)+random(4);
    end
    
    num_FAP = input('How many FAPS?');
    dx_FAP = zeros(1,num_FAP);
    dy_FAP = zeros(1,num_FAP);
    num_femto = 2 * num_FAP;
    dx_femto = zeros(1,num_femto);
    dy_femto = zeros(1,num_femto);
    
    for i=1:num_FAP
        dx_FAP(i) = x(i); % the matrix with teh FAPs' coordinates
        dy_FAP(i) = y(i);
        dx_femto(2*i-1) = xf(2*i-1);% NRT users
        dy_femto(2*i-1) = yf(2*i-1);
        dx_femto(2*i) = xf(2*i);% RT users
        dy_femto(2*i) = yf(2*i);
    end
    
    num_macro = input('How many macro users?');
    num_nrt = input('How many nrt macro users?');
    %num_rt = num_macro-num_nrt;
    dx_macro = randi([100 900],1,num_macro);
    dy_macro = randi([100 900],1,num_macro);
    
    % Matrix of all distances between stations and path gains
    num_users = num_macro + num_femto;  
    
    D = zeros(num_FAP+1,num_users); % the +1 is to count the MBS (the Macrocell base station) 
    % rows = base stations, columns = users
    Gain = zeros(num_FAP+1,num_users);
    kc = 0.1;% fixed propagation loss during cellular transmissions to MBS
    kf = 0.01;% fixed loss between femto-user i to their FAP(i).
    w = 0.3162;% partition loss during indoor-to-outdoor propagation
    w2 = 0.1;% w^2 (double loss)
    a = 4;% outdoor path loss exponent
    b = 3;% indoor path loss exponent
    
    for i=1:num_FAP+1
        for j=1:num_users
            if (i==1)&&(j<=num_macro) % m -> MBS
                D(i,j) = euklidian(dx_macro(j),dy_macro(j),500,500);
                Gain(i,j) = kc * power(D(i,j),-a);
            end
            if (i>1)&&(j<=num_macro) % m -> FAPs
                ifap = i-1;
                D(i,j) = euklidian(dx_macro(j),dy_macro(j),dx_FAP(ifap),dy_FAP(ifap));
                Gain(i,j) = kc * w * power(D(i,j),-a);
            end
            if (i==1)&&(j>num_macro) % f -> MBS
                k = j - num_macro;
                D(i,j) = euklidian(dx_femto(k),dy_femto(k),500,500);
                Gain(i,j) = kc * w * power(D(i,j),-a);
            end
            if (i>1)&&(j>num_macro) % f -> FAPs
                k = j - num_macro; % the femto-user
                ifap = i-1; % which FAP we check 
                D(i,j) = euklidian(dx_femto(k),dy_femto(k),dx_FAP(ifap),dy_FAP(ifap)); 
                if (2*ifap-1 == k)||(2*ifap == k) % f -> FAP f belongs to that FAP - he is a home user!
                    Gain(i,j) = kf * power(D(i,j),-b);
                else % f user is from another FAP (a neighboor!)
                    Gain(i,j) = kc * w2 * power(D(i,j),-a);          
                end
            end
        end
    end
    
    % Initializations
    Pmax = 2.0; % in Watt
    Rmax = 2.4*10^6; % in bps
    Rtarget = 64*10^3; % in bps
    MF = 10*10^3; % in bps
    Rmaxi = Rtarget + MF;
    Rmini = Rtarget - MF;
    W = 10^6; %bandwidth in Hz
    noise = 5 * 10^(-16); %noise
    e = 10^(-5);
    G = W/Rmaxi;
    
    % Power Control Algorithm

    c = 10^10; %femto pricing factor
    p = 0:0.00001:2;
    p_users = 2*ones(1,num_users);
    % in order of: nrt_m,rt_m,2 femto_users(nrt,rt)/FAP
    converged = false;
    k = 1;
    new_p_users = p_users;
    isxis = cell(1,num_users);
    
    while (converged == false) && (k < 80)
        I = Gain*p_users'; % total intereference in every BS (num_FAP+1 X 1)
        for i=1:num_users
            if i <= num_macro % MACRO USERS
               I_i = I(1) - Gain(1,i) * p_users(i) + noise;
               if i <= num_nrt % NRT macro users
                  U = ((log(0.001*Rmax*(power((1-exp(-3.7*(G*(Gain(1,i)*p)/I_i))),80))+1))./p);
                  %find the max of U
                  Um = U(~isnan(U)); 
                  if (1 && all(Um == 0)) 
                     root = Pmax;
                  else
                     [~,in] = max(U);
                     root = p(in);
                  end
               else % RT macro users
                  U = (power(1-exp(-(Rmaxi*(power((1-exp(-3.7*(G*(Gain(1,i)*p)/I_i))),80))-Rmini)),2000)./p);
                  %find the max of U 
                  l1 = length(U);
                  Um = U(~isinf(U));
                  l2 = length(Um);
                  if (l2 == 0)
                     root = Pmax;
                  else
                     [~,in] = max(U(~isinf(U)));
                     ini = l1-l2+in;
                     root = p(ini);
                  end
               end
            else % FEMTO USERS (i > num_macro)
               j = i - num_macro;
               if mod(j,2) == 1 % NRT femto users
                  f = (j+1)/2; % f is the FAP that the j-th femto user belongs
                  I_i = I(f+1) - Gain(f+1,i) * p_users(i) + noise;
                  U = ((log(0.001*Rmax*(power((1-exp(-3.7*(G*(Gain(f+1,i)*p)/I_i))),80))+1))./p); %-c*(exp(p)-1);
                  %find the max of U
                  Um = U(~isnan(U)); 
                  if (1 && all(Um == 0)) 
                     root = Pmax;
                  else
                     [~,in] = max(U);
                     root = p(in);
                  end
               else % RT femto users
                  f = j/2;
                  I_i = I(f+1) - Gain(f+1,i) * p_users(i) + noise;
                  U = (power(1-exp(-(Rmaxi*(power((1-exp(-3.7*(G*(Gain(f+1,i)*p)/I_i))),80))-Rmini)),2000)./p)-c*(exp(p)-1);
                  %find the max of U
                  l1 = length(U);
                  Um = U(~isinf(U));
                  l2 = length(Um);
                  if (l2 == 0)
                     root = Pmax;
                  else 
                     [~,in] = max(U(~isinf(U)));
                     ini = l1-l2+in;
                     root = p(ini);
                  end
               end    
            end   
            new_p_users(i) = root;
        end
        
        % Test Convergence
        count = 0;
        for i=1:num_users
           if (abs(new_p_users(i) - p_users(i)) < e)
              count = count + 1;
           end
        end
        if count == num_users
           converged = true;
        else
           converged = false;
        end
        for i = 1:num_users
          isxis{i}(k) = p_users(i);       
        end 
        p_users = new_p_users
        k = k+1
    end
    
    % Calculate the Users' Utilities 
    
    I = Gain*p_users';
    fgama = zeros(1,num_users);
    Rate = zeros(1,num_users);
    Util = zeros(1,num_users);
    
    for i=1:num_users
       if i <= num_macro % MACRO USERS
          I_i = I(1) - Gain(1,i) * p_users(i) + noise;
          if i <= num_nrt % NRT macro users
             fgama(i) = power((1-exp(-3.7*(G*(Gain(1,i)*p_users(i))/I_i))),80);
             Rate(i) = fgama(i) * Rmax;
             Util(i) = (log(0.001*Rmax*fgama(i)+1))./p_users(i);
          else % RT macro users
             fgama(i) = power((1-exp(-3.7*(G*(Gain(1,i)*p_users(i))/I_i))),80);
             Rate(i) = fgama(i) * Rmaxi;
             Util(i) = power(1-exp(-(Rmaxi*fgama(i)-Rmini)),2000)./p_users(i);      
          end
       else % FEMTO USERS (i > num_macro)
          j = i - num_macro;
          if mod(j,2) == 1 % NRT femto users  
             f = (j+1)/2; 
             I_i = I(f+1) - Gain(f+1,i) * p_users(i) + noise;
             fgama(i) = power((1-exp(-3.7*(G*(Gain(f+1,i)*p_users(i))/I_i))),80);
             Rate(i) = fgama(i) * Rmax;
             Util(i) = ((log(0.001*Rmax*fgama(i)+1))./p_users(i));
          else % RT femto users
             f = j/2;
             I_i = I(f+1) - Gain(f+1,i) * p_users(i) + noise;
             fgama(i) = power((1-exp(-3.7*(G*(Gain(f+1,i)*p_users(i))/I_i))),80);
             Rate(i) = fgama(i) * Rmaxi;
             Util(i) = (power(1-exp(-(Rmaxi*fgama(i)-Rmini)),2000)./p_users(i));
          end           
       end
    end
    
    figure;
    % Power vs iterations for macro users
    iterations = 1:k-1;
    for i = 1:num_macro
          hold on;
          title('P vs iterations for macro users');
          if i <= num_nrt
             plot(iterations,isxis{i},'r*');
          else
             plot(iterations,isxis{i},'b*');
          end
    end
    xlabel('Iterations');
    ylabel('User Transmission Power in Watt');
    hold off;
    
    figure;
    % Power vs iterations for femto users
    iterations = 1:k-1;
    for i = num_macro+1:num_users
          hold on;
          title('P vs iterations for femto users');
          j = i - num_macro;
          if mod(j,2) == 1
             plot(iterations,isxis{i},'r*');
          else
             plot(iterations,isxis{i},'b*');
          end
    end
    xlabel('Iterations');
    ylabel('User Transmission Power in Watt');
    hold off;
    
    id = 1:num_users;
    figure;bar(id,fgama,'g');title('f(γ)');
    figure;bar(id,Rate,'r');title('Users Rate');
    figure;bar(id,Util);title('Users Utilities');
  
    % Draw the Macrocell with all the users
    figure;
    draw(dx_macro,dy_macro,dx_femto,dy_femto,dx_FAP,dy_FAP,num_nrt,num_macro,num_femto);
end

function d = euklidian(x,y,z,w)
    d = sqrt((x-z).^2+(y-w).^2);
end

function draw(dx_macro,dy_macro,dx_femto,dy_femto,dx_FAP,dy_FAP,num_nrt,num_macro,num_femto)
    
   mbsx = 500; % MBS coordinates
   mbsy = 500;
   num_rt = num_macro - num_nrt;

   % distinguish macro users to RT and NRT
   dx_macro_nrt = zeros(1,num_nrt);
   dy_macro_nrt = zeros(1,num_nrt);
   dx_macro_rt = zeros(1,num_rt);
   dy_macro_rt = zeros(1,num_rt);
   for i=1:num_macro
       if i<=num_nrt
           dx_macro_nrt(i) = dx_macro(i);
           dy_macro_nrt(i) = dy_macro(i);
       else
           dx_macro_rt(i) = dx_macro(i);
           dy_macro_rt(i) = dy_macro(i);           
       end
   end
   
   % distinguish femto users to RT and NRT
   num_FAP = num_femto/2;
   dx_femto_nrt = zeros(1,num_FAP);
   dy_femto_nrt = zeros(1,num_FAP);
   dx_femto_rt = zeros(1,num_FAP);
   dy_femto_rt = zeros(1,num_FAP);
   
   k = 1;
   w = 1;
   for i=1:num_femto
      if mod(i,2) == 1 % 2 users/FAP, first is a NRT one, second is a RT
          dx_femto_nrt(k) = dx_femto(i);
          dy_femto_nrt(k) = dy_femto(i);
          k = k+1;
      else
          dx_femto_rt(w) = dx_femto(i);
          dy_femto_rt(w) = dy_femto(i);
          w = w+1;
      end
   end
   
   hold on; % red = NRT, blue = RT
   grid on;
   plot(mbsx,mbsy,'kd',dx_FAP,dy_FAP,'bd',dx_macro_nrt,dy_macro_nrt,'rx',dx_macro_rt,dy_macro_rt,'bx',dx_femto_nrt,dy_femto_nrt,'ro',dx_femto_rt,dy_femto_rt,'bo');
   xlim([0 1000]);ylim([0 1000]);
   circle(mbsx,mbsy,500,'r');
   for i=1:num_FAP
      circle(dx_FAP(i),dy_FAP(i),50,'g');
   end
   axis equal;
   hold off;
end

function circle(x,y,r,color)
   th = 0:pi/500:2*pi;
   xunit = r * cos(th) + x;
   yunit = r * sin(th) + y;
   plot(xunit, yunit, color);
end
