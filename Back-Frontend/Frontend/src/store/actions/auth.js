import http from '../../services/httpService.js';
import {
  LOGIN_REQUEST,
  LOGIN_SUCCESS,
  LOGIN_FAILURE,
  LOGOUT_REQUEST,
} from '../constants/auth';

const apiEndpoint = 'https://reqres.in';

export const loginRequest = () => ({
  type: LOGIN_REQUEST,
});

export const loginSuccess = ({ token }) => {
  console.log(token);

  localStorage.setItem('token', token);
  return {
    type: LOGIN_SUCCESS,
    token: token,
    username: 'ketan',
  };
};

export const loginFailure = (error) => {
  console.log(error);

  return {
    type: LOGIN_FAILURE,
    status: 'error.response.status',
    statusText: 'login failed',
  };
};

export const logout = () => {
  localStorage.removeItem('token');
  return { type: LOGOUT_REQUEST };
};

export const loginUser = (email, password) => async (dispatch) => {
  console.log({ email, password });

  dispatch(loginRequest());
  try {
    const { data } = await http.post(`https://reqres.in/api/login`, {
      email,
      password,
    });
    //decode token then
    dispatch(loginSuccess(data));
  } catch (error) {
    dispatch(loginFailure(error));
  }
};
