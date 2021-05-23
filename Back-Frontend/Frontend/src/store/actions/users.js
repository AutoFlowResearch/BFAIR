import {
  RECEIVE_USERS,
  REQUEST_USERS,
  FAILURE_USERS,
  RECEIVE_USER,
  REQUEST_USER,
  FAILURE_USER,
  UPDATE_USER_DATA,
  CLEAR_USER_DATA,
} from '../constants/users';
import { userslist, getUserById as userById } from '../../services/users';

// list of users
const requestUsers = () => ({
  type: REQUEST_USERS,
});

const receiveUsers = (data) => ({
  type: RECEIVE_USERS,
  data,
});

const failureUsers = (error) => ({
  type: FAILURE_USERS,
  error,
});

// user by id

const requestUser = () => ({
  type: REQUEST_USER,
});

const receiveUser = (data) => ({
  type: RECEIVE_USER,
  data,
});

const failureUser = (error) => ({
  type: FAILURE_USER,
  error,
});

export const updateData = (payload) => ({
  type: UPDATE_USER_DATA,
  payload,
});

export const clearUserData = () => ({
  type: CLEAR_USER_DATA,
});

export const getUsers = () => async (dispatch) => {
  dispatch(requestUsers());
  try {
    const data = await userslist();
    dispatch(receiveUsers(data));
  } catch (error) {
    dispatch(failureUsers(error));
  }
};

export const getUserByid = (id) => async (dispatch) => {
  dispatch(requestUser());
  try {
    const data = await userById(id);
    dispatch(receiveUser(data));
  } catch (error) {
    dispatch(failureUser(error));
  }
};
