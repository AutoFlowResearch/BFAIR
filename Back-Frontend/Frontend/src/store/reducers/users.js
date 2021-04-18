import {
  RECEIVE_USERS,
  REQUEST_USERS,
  FAILURE_USERS,
  RECEIVE_USER,
  REQUEST_USER,
  FAILURE_USER,
  UPDATE_USER_DATA,
  CLEAR_USER_DATA
} from "../constants/users";

const initialState = {
  userslist: {
    data: [],
    isFetching: false,
    error: null
  },

  updateuser: {
    isFetching: false,
    error: null,
    data: {
      id: "",
      name: "",
      username: "",
      email: "",
      address: {
        street: "",
        suite: "",
        city: "",
        zipcode: "",
        geo: {
          lat: "",
          lng: ""
        }
      },
      phone: "",
      website: "",
      company: {
        name: "",
        catchPhrase: "",
        bs: ""
      }
    }
  }
};

const userslist = (state, action) => {
  switch (action.type) {
    case REQUEST_USERS:
      return { ...state, isFetching: true };
    case RECEIVE_USERS:
      return { ...state, isFetching: false, data: action.data };
    case FAILURE_USERS:
      return { ...state, error: action.error };
    default:
      return state;
  }
};

const updateuser = (state, action) => {
  switch (action.type) {
    case REQUEST_USER:
      return { ...state, isFetching: true };
    case RECEIVE_USER:
      return { ...state, data: action.data, isFetching: false };
    case FAILURE_USER:
      return { ...state, error: action.error };
    case UPDATE_USER_DATA:
      return { ...state, data: action.payload };
    case CLEAR_USER_DATA:
      return initialState.updateuser;
    default:
      return state;
  }
};

export const usermanagement = (state = initialState, action) => {
  // console.log(state);
  switch (action.type) {
    default:
      return {
        userslist: userslist(state.userslist, action),
        updateuser: updateuser(state.updateuser, action)
      };
  }
};
