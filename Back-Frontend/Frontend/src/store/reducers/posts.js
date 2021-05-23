import lodash from 'lodash';

import {
  RECEIVE_POSTS,
  REQUEST_POSTS,
  FAILURE_POSTS,
  RECEIVE_POST,
  REQUEST_POST,
  FAILURE_POST,
  UPDATE_POST_DATA,
  CLEAR_POST_DATA,
  CREATE_POST,
} from '../constants/posts';

const initialState = {
  postsList: {
    data: [],
    isFetching: false,
    error: null,
  },

  updatePost: {
    isFetching: false,
    error: null,
    data: {
      title: '',
      body: '',
    },
  },
};

const postsList = (state, action) => {
  switch (action.type) {
    case REQUEST_POSTS:
      return { ...state, isFetching: true };
    case RECEIVE_POSTS:
      return {
        ...state,
        isFetching: false,
        data: lodash.unionBy(state.data, action.data, 'id'),
      };
    case FAILURE_POSTS:
      return { ...state, error: action.error };
    case CREATE_POST:
      const data = [...state.data, action.payload];
      return { ...state, data };
    default:
      return state;
  }
};

const updatePost = (state, action) => {
  switch (action.type) {
    case REQUEST_POST:
      return { ...state, isFetching: true };
    case RECEIVE_POST:
      return { ...state, data: action.data, isFetching: false };
    case FAILURE_POST:
      return { ...state, error: action.error };
    case UPDATE_POST_DATA:
      return { ...state, data: action.payload };
    case CLEAR_POST_DATA:
      return initialState.updatePost;
    default:
      return state;
  }
};

export const postManagement = (state = initialState, action) => {
  // console.log(state);
  switch (action.type) {
    default:
      return {
        postsList: postsList(state.postsList, action),
        updatePost: updatePost(state.updatePost, action),
      };
  }
};
