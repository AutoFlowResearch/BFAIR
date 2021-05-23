import {
  LOGIN_REQUEST,
  LOGIN_SUCCESS,
  LOGIN_FAILURE,
  LOGOUT_REQUEST,
} from '../constants/auth';

const isuser = localStorage.getItem('token');

const innitialState = {
  token: isuser || null,
  userName: isuser ? 'lawrence' : null,
  isAuthenticated: !!isuser,
  isAuthenticating: false,
  statusText: isuser ? 'Logged In' : null,
};

export const auth = (state = innitialState, action) => {
  switch (action.type) {
    case LOGIN_REQUEST:
      return { ...state, isAuthenticating: true };
    case LOGIN_SUCCESS:
      return {
        ...state,
        token: action.token,
        userName: action.userName,
        isAuthenticating: false,
        isAuthenticated: true,
        statusText: 'You have been successfully logged in',
      };
    case LOGIN_FAILURE:
      return {
        ...state,
        token: null,
        userName: null,
        isAuthenticated: false,
        isAuthenticating: false,
        statusText: `Authentication Error ${action.status} ${action.statusText}`,
      };
    case LOGOUT_REQUEST:
      return {
        ...state,
        isAuthenticated: false,
        token: null,
        userName: null,
        statusText: 'You have been successfully logged out.',
      };
    default:
      return state;
  }
};
